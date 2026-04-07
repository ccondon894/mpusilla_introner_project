#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import logging
import os

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'results', 'genotyping', 'logs', 'introner_detection.log'))
    ]
)
logger = logging.getLogger(__name__)

def load_input_files(bed_file, depth_file, call_file):
    """Load and prepare input files."""
    logger.info(f"Loading input files...")
    
    # Load depth data
    depth_df = pd.read_csv(depth_file, sep="\t", header=None, 
                          names=["chrom", "pos", "depth"])
    depth_df = depth_df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    
    # Load call data
    gt_matrix_df = pd.read_csv(call_file, sep="\t", header=0)
    
    # Load bed data - detect if first line is a header
    with open(bed_file, 'r') as f:
        first_line = f.readline().strip()
        # Check if first line looks like a header (contains 'chrom', 'start', or 'ortholog')
        if any(x in first_line.lower() for x in ['chrom', 'start', 'ortholog']):
            logger.info("BED file appears to have a header, reading accordingly")
            bed_df = pd.read_csv(bed_file, sep="\t")
        else:
            # No header, read as usual
            bed_df = pd.read_csv(bed_file, sep="\t", header=None, 
                                names=["chrom", "start", "end", "sample", "ortholog_id",
                                       "sequence_id", "family", "gene", "splice_site", "presence", "orientation",
                                       "left_reverse", "right_reverse"])
    
    return bed_df, depth_df, gt_matrix_df

def get_region_coverage(depth_df, chrom, start, end):
    """Calculate coverage statistics for a genomic region."""
    # Filter depth data for this region
    region_cov = depth_df[(depth_df["chrom"] == chrom) & 
                          (depth_df["pos"] >= start) & 
                          (depth_df["pos"] < end)]
    
    # Get all positions in the range
    all_positions = set(range(start, end))
    covered_positions = set(region_cov["pos"].unique())
    missing_positions = all_positions - covered_positions
    
    # Calculate coverage metrics
    total_bases = end - start
    covered_bases = len(covered_positions)
    coverage_ratio = covered_bases / total_bases if total_bases > 0 else 0
    
    # Calculate depth metrics
    sum_depth = region_cov["depth"].sum()
    mean_depth = sum_depth / total_bases if total_bases > 0 else 0
    
    # Calculate median depth (more robust to outliers)
    if not region_cov.empty:
        median_depth = region_cov["depth"].median()
    else:
        median_depth = 0
    
    # Count positions with zero depth
    zero_depth_count = total_bases - covered_bases
    
    # Calculate uniformity measure
    if not region_cov.empty and mean_depth > 0:
        uniformity = median_depth / mean_depth
    else:
        uniformity = 0
    
    return {
        "total_bases": total_bases,
        "covered_bases": covered_bases,
        "coverage_ratio": coverage_ratio,
        "mean_depth": mean_depth,
        "median_depth": median_depth,
        "zero_depth_count": zero_depth_count,
        "uniformity": uniformity
    }

def analyze_spanning_reads(depth_df, chrom, junction_pos, window=10):
    """
    Analyze coverage specifically at the junction between flanking and introner regions.
    Looking for evidence of reads that span both regions.
    """
    left_side = junction_pos - window
    right_side = junction_pos + window
    
    # Get coverage around the junction
    junction_cov = depth_df[(depth_df["chrom"] == chrom) & 
                           (depth_df["pos"] >= left_side) & 
                           (depth_df["pos"] < right_side)]
    
    # Check for continuous coverage across the junction
    all_positions = set(range(left_side, right_side))
    covered_positions = set(junction_cov["pos"].unique())
    missing_positions = all_positions - covered_positions
    
    has_continuous_coverage = len(missing_positions) == 0
    
    # Calculate depth change across junction
    left_depth = junction_cov[junction_cov["pos"] < junction_pos]["depth"].mean() if not junction_cov.empty else 0
    right_depth = junction_cov[junction_cov["pos"] >= junction_pos]["depth"].mean() if not junction_cov.empty else 0
    depth_ratio = right_depth / left_depth if left_depth > 0 else 0
    
    return {
        "has_continuous_coverage": has_continuous_coverage,
        "left_depth": left_depth,
        "right_depth": right_depth,
        "depth_ratio": depth_ratio,
        "missing_positions": list(missing_positions)
    }

def model_coverage_pattern(depth_df, chrom, start, end):
    """
    Model the coverage pattern to determine if it matches expected introner patterns.
    """
    region_cov = depth_df[(depth_df["chrom"] == chrom) & 
                          (depth_df["pos"] >= start) & 
                          (depth_df["pos"] < end)]
    
    # If no coverage data, return early
    if region_cov.empty:
        return {
            "best_pattern": "no_data",
            "pattern_scores": {"no_data": 1},
            "middle_flank_ratio": 0,
            "cv": 0
        }
    
    # Divide region into bins for pattern analysis
    bins = 10
    bin_size = max(1, (end - start) // bins)
    bin_depths = []
    
    for i in range(bins):
        bin_start = start + i * bin_size
        bin_end = min(bin_start + bin_size, end)
        bin_data = region_cov[(region_cov["pos"] >= bin_start) & 
                            (region_cov["pos"] < bin_end)]
        bin_depths.append(bin_data["depth"].mean() if not bin_data.empty else 0)
    
    # Calculate pattern metrics
    mean_depth = np.mean(bin_depths) if bin_depths else 0
    std_depth = np.std(bin_depths) if bin_depths else 0
    cv = std_depth / mean_depth if mean_depth > 0 else 0  # Coefficient of variation
    
    # Middle vs flanks ratio
    middle_idx = bins // 3
    flank_idx = 2 * bins // 3
    
    left_flank = np.mean(bin_depths[:middle_idx]) if bin_depths else 0
    middle = np.mean(bin_depths[middle_idx:flank_idx]) if bin_depths else 0
    right_flank = np.mean(bin_depths[flank_idx:]) if bin_depths else 0
    
    middle_flank_ratio = middle / ((left_flank + right_flank) / 2) if (left_flank + right_flank) > 0 else 0
    
    # Pattern matching
    pattern_scores = {
        "introner": 0,
        "continuous": 0,
        "junction": 0,
        "absent": 0
    }
    
    # Higher middle:flank ratio favors introner pattern (potential multimap)
    if middle_flank_ratio > 1.5:
        pattern_scores["introner"] += 1
    
    # Low CV favors continuous pattern
    if cv < 0.2:
        pattern_scores["continuous"] += 1
    
    # Check for sharp transitions at junctions
    if middle_idx > 0 and flank_idx > 0 and len(bin_depths) > flank_idx:
        left_junction = abs(bin_depths[middle_idx] - bin_depths[middle_idx-1]) / mean_depth if mean_depth > 0 else 0
        right_junction = abs(bin_depths[flank_idx] - bin_depths[flank_idx-1]) / mean_depth if mean_depth > 0 else 0
        
        if left_junction > 0.5 or right_junction > 0.5:
            pattern_scores["junction"] += 1
    
    # Very low coverage suggests absence
    if mean_depth < 5:
        pattern_scores["absent"] += 1
    
    # Get best matching pattern
    best_pattern = max(pattern_scores.items(), key=lambda x: x[1])[0]
    
    return {
        "best_pattern": best_pattern,
        "pattern_scores": pattern_scores,
        "middle_flank_ratio": middle_flank_ratio,
        "cv": cv
    }

def call_presence(row, call, depth_df):
    """Improved introner detection focused on spanning reads."""
    chrom = row["contig"]
    start = int(row["start"])
    end = int(row["end"])
    
    # Define regions
    left_start, left_end = start, start + 100
    right_start, right_end = end - 100, end
    middle_start, middle_end = left_end, right_start
    
    # Get coverage for flanking regions
    left_cov = get_region_coverage(depth_df, chrom, left_start, left_end)
    right_cov = get_region_coverage(depth_df, chrom, right_start, right_end)
    
    # Check if flanking regions have sufficient coverage (minimum mapping threshold)
    if left_cov["mean_depth"] < 10 or right_cov["mean_depth"] < 10:
        logger.debug(f"Insufficient flanking coverage: {left_cov['mean_depth']}, {right_cov['mean_depth']}")
        return left_cov["mean_depth"], 0, right_cov["mean_depth"], 3
    
   
    # For presence calls, we need to validate with evidence
    middle_cov = get_region_coverage(depth_df, chrom, middle_start, middle_end)
    
    # Ensure adequate coverage across introner region
    min_acceptable_coverage = 0.8  # At least 80% of positions must have coverage
    if middle_cov["coverage_ratio"] < min_acceptable_coverage:
        logger.debug(f"Insufficient middle coverage ratio: {middle_cov['coverage_ratio']}")
        return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], 2
    
    # Check for minimum depth
    if middle_cov["mean_depth"] < 10:
        logger.debug(f"Low middle mean depth: {middle_cov['mean_depth']}")
        return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], 2
    
    # Most importantly, check for reads spanning the junctions
    left_junction = analyze_spanning_reads(depth_df, chrom, middle_start)
    right_junction = analyze_spanning_reads(depth_df, chrom, middle_end)
    
    # Critical test: do we have evidence of continuous coverage across junctions?
    if not left_junction["has_continuous_coverage"] or not right_junction["has_continuous_coverage"]:
        logger.debug(f"Missing spanning reads: left={left_junction['has_continuous_coverage']}, right={right_junction['has_continuous_coverage']}")
        return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], 2
    
    # Check for abnormal depth differences that might indicate multimap issues
    if middle_cov["mean_depth"] > 3 * max(left_cov["mean_depth"], right_cov["mean_depth"]):
        # Suspiciously high middle coverage might indicate multimapping
        logger.debug(f"Suspiciously high middle coverage: middle={middle_cov['mean_depth']}, flanks={max(left_cov['mean_depth'], right_cov['mean_depth'])}")
        return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], 2
    
    # Check coverage uniformity across the introner
    if middle_cov["uniformity"] < 0.6:  # Very uneven coverage
        logger.debug(f"Uneven middle coverage: uniformity={middle_cov['uniformity']}")
        return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], 2
    
    # Model the overall coverage pattern
    pattern = model_coverage_pattern(depth_df, chrom, start, end)
    if pattern["best_pattern"] == "absent":
        logger.debug(f"Coverage pattern suggests absence: {pattern}")
        return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], 2
    
    # Passed all tests - valid introner presence
    return left_cov["mean_depth"], middle_cov["mean_depth"], right_cov["mean_depth"], call

def main():
    # Parse command line arguments
    if len(sys.argv) != 5:
        print("Usage: python call_locus_from_bam_coverage.py bed_file depth_file call_file output_bed_file")
        sys.exit(1)
    
    bed_file = sys.argv[1]
    depth_file = sys.argv[2]
    gt_matrix_file = sys.argv[3]
    output_bed_file = sys.argv[4]

    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_bed_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load input files
    bed_df, depth_df, gt_df = load_input_files(bed_file, depth_file, gt_matrix_file)
    
    # Process each locus
    total_loci = len(bed_df)
    changed_calls = 0
    
    logger.info(f"Processing {total_loci} loci...")
    
    with open(output_bed_file, 'w') as o:
        # Write header
        header = "contig\tstart\tend\tsample\tortholog_id\tsequence_id\tfamily\tgene\tsplice_site\torientation\tleft_reverse\tright_reverse\tleft_mean_depth\tmean_middle_depth\tright_mean_depth\tcall\tnew_call\n"
        o.write(header)
        
        # Process each locus
        for index, row in bed_df.iterrows():

            ortholog_id = row["ortholog_id"]
            sample_id = row["sample"]
            print(ortholog_id, sample_id)
            # Log progress every 100 loci
            if index % 100 == 0:
                logger.info(f"Processing locus {index+1}/{total_loci}: {ortholog_id}, {sample_id}")
            
            # Get original call
            call = gt_df.loc[(gt_df['ortholog_id'] == ortholog_id) & (gt_df['sample'] == sample_id), 'presence'].iloc[0]
            
            # Call presence/absence
            left_mean_depth, mean_middle_depth, right_mean_depth, new_call = call_presence(row, call, depth_df)
            
            # Write output
            data_row = f"{row['contig']}\t{row['start']}\t{row['end']}\t{row['sample']}\t{row['ortholog_id']}" \
                     + f"\t{row['sequence_id']}\t{row['family']}\t{row['gene']}\t{row['splice_site']}" \
                     + f"\t{row['orientation']}\t{row['left_reverse']}\t{row['right_reverse']}" \
                     + f"\t{left_mean_depth}\t{mean_middle_depth}\t{right_mean_depth}\t{call}\t{new_call}\n"
            o.write(data_row)
            
            # Log changes
            if call != new_call:
                changed_calls += 1
                logger.info(f"Changed call for {ortholog_id}, {sample_id}: {call} -> {new_call}")
                # update gt matrix
                

    logger.info(f"Processing complete. Changed {changed_calls}/{total_loci} calls.")
    logger.info(f"Results written to {output_bed_file}")

if __name__ == "__main__":
    main()