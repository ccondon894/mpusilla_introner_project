#!/usr/bin/env python3
import argparse
import gffutils
import pandas as pd
import numpy as np
from collections import defaultdict
import logging
from pathlib import Path
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_target_id(attributes):
    """Extract target ID from GTF attributes"""
    # Find the Target attribute and extract the first part (ID)
    target_match = re.search(r'Target "([^"]+)"', attributes)
    if target_match:
        # Extract just the ID part (before the coordinates)
        return target_match.group(1).split()[0]
    return None

class GTFStats:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.db = self._create_db()
        
    def _create_db(self):
        """Create a GFFUtils database from GTF file"""
        def custom_id_spec(feature):
            """Custom ID specification to handle malformed ID fields"""
            # If feature has an ID attribute, clean it up
            if 'ID' in feature.attributes:
                # Take the first ID value and clean it
                raw_id = feature.attributes['ID'][0]
                # Remove trailing commas and spaces
                clean_id = raw_id.rstrip(', ')
                return clean_id

            # Fall back to default behavior
            return feature.id

        return gffutils.create_db(
            self.gtf_file,
            dbfn=':memory:',
            force=True,
            merge_strategy='create_unique',
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            id_spec=custom_id_spec
        )
    
    def get_basic_stats(self):
        """Get basic statistics about the GTF"""
        stats = {
            'gene_count': len(list(self.db.features_of_type('gene'))),
            'transcript_count': len(list(self.db.features_of_type('transcript'))),
            'exon_count': len(list(self.db.features_of_type('exon'))),
            'cds_count': len(list(self.db.features_of_type('CDS')))
        }
        
        # Calculate gene length statistics
        gene_lengths = []
        gene_identities = []  # Add collection of identity scores
        for gene in self.db.features_of_type('gene'):
            gene_lengths.append(gene.end - gene.start + 1)
            
            # Extract identity score if available
            identity_str = gene.attributes.get('Identity', ['0'])[0].strip('"')
            try:
                identity = float(identity_str)
                gene_identities.append(identity)
            except ValueError:
                pass
            
        if gene_lengths:
            stats.update({
                'mean_gene_length': np.mean(gene_lengths),
                'median_gene_length': np.median(gene_lengths),
                'min_gene_length': np.min(gene_lengths),
                'max_gene_length': np.max(gene_lengths)
            })
            
        if gene_identities:
            stats.update({
                'mean_identity': np.mean(gene_identities),
                'median_identity': np.median(gene_identities),
                'min_identity': np.min(gene_identities),
                'max_identity': np.max(gene_identities)
            })
        
        # Calculate exons per gene
        exons_per_gene = defaultdict(int)
        for exon in self.db.features_of_type('exon'):
            parent_gene = list(self.db.parents(exon, featuretype='gene'))
            if parent_gene:
                exons_per_gene[parent_gene[0].id] += 1
                
        if exons_per_gene:
            stats.update({
                'mean_exons_per_gene': np.mean(list(exons_per_gene.values())),
                'median_exons_per_gene': np.median(list(exons_per_gene.values())),
                'min_exons_per_gene': np.min(list(exons_per_gene.values())),
                'max_exons_per_gene': np.max(list(exons_per_gene.values()))
            })
            
        return stats
    
    def get_chromosome_distribution(self):
        """Get distribution of features across chromosomes/scaffolds"""
        chrom_dist = defaultdict(lambda: defaultdict(int))
        
        for feat_type in ['gene', 'transcript', 'exon', 'CDS']:
            for feature in self.db.features_of_type(feat_type):
                chrom_dist[feature.seqid][feat_type] += 1
                
        return dict(chrom_dist)
    
    def get_target_mapping(self):
        """Get mapping of genes to their target IDs"""
        target_map = {}
        for gene in self.db.features_of_type('gene'):
            target_id = parse_target_id(str(gene.attributes))
            if target_id:
                target_map[gene.id] = target_id
        return target_map

def compare_gtfs(ref_gtf, query_gtf):
    """Compare two GTF files and generate comparison statistics"""
    logger.info(f"Analyzing reference GTF: {ref_gtf}")
    ref_stats = GTFStats(ref_gtf)
    ref_basic = ref_stats.get_basic_stats()
    ref_chrom = ref_stats.get_chromosome_distribution()
    ref_targets = set(ref_stats.get_target_mapping().values())  # Get reference target IDs
    
    logger.info(f"Analyzing query GTF: {query_gtf}")
    query_stats = GTFStats(query_gtf)
    query_basic = query_stats.get_basic_stats()
    query_chrom = query_stats.get_chromosome_distribution()
    query_targets = set(query_stats.get_target_mapping().values())  # Get query target IDs
    
    # Calculate overlap statistics using target IDs
    shared_targets = ref_targets.intersection(query_targets)
    unique_to_ref = ref_targets - query_targets
    unique_to_query = query_targets - ref_targets
    
    overlap_stats = {
        'shared_genes': len(shared_targets),
        'unique_to_ref': len(unique_to_ref),
        'unique_to_query': len(unique_to_query),
        'percent_ref_recovered': len(shared_targets) / len(ref_targets) * 100 if ref_targets else 0
    }
    
    # Get identity score distribution for shared genes
    shared_gene_stats = {}
    if shared_targets:
        identity_scores = []
        for gene in query_stats.db.features_of_type('gene'):
            target_id = parse_target_id(str(gene.attributes))
            if target_id in shared_targets:
                identity_str = gene.attributes.get('Identity', ['0'])[0].strip('"')
                try:
                    identity = float(identity_str)
                    identity_scores.append(identity)
                except ValueError:
                    continue
        
        if identity_scores:
            shared_gene_stats.update({
                'mean_shared_identity': np.mean(identity_scores),
                'median_shared_identity': np.median(identity_scores),
                'min_shared_identity': np.min(identity_scores),
                'max_shared_identity': np.max(identity_scores)
            })
    
    return {
        'reference_stats': ref_basic,
        'query_stats': query_basic,
        'reference_chromosome_dist': ref_chrom,
        'query_chromosome_dist': query_chrom,
        'overlap_statistics': overlap_stats,
        'shared_gene_stats': shared_gene_stats
    }

def write_stats(stats, output_prefix):
    """Write statistics to files"""
    # Basic stats comparison
    # Ensure both dictionaries have the same keys
    all_keys = set(stats['reference_stats'].keys()) | set(stats['query_stats'].keys())
    ref_values = [stats['reference_stats'].get(key, 0) for key in all_keys]
    query_values = [stats['query_stats'].get(key, 0) for key in all_keys]

    basic_comp = pd.DataFrame({
        'Metric': list(all_keys),
        'Reference': ref_values,
        'Query': query_values
    })
    basic_comp.to_csv(f"{output_prefix}_basic_stats.tsv", sep='\t', index=False)
    
    # Chromosome distribution
    ref_chrom = pd.DataFrame(stats['reference_chromosome_dist']).fillna(0)
    query_chrom = pd.DataFrame(stats['query_chromosome_dist']).fillna(0)
    
    ref_chrom.to_csv(f"{output_prefix}_ref_chrom_dist.tsv", sep='\t')
    query_chrom.to_csv(f"{output_prefix}_query_chrom_dist.tsv", sep='\t')
    
    # Overlap statistics
    overlap = pd.DataFrame({
        'Metric': list(stats['overlap_statistics'].keys()) + list(stats['shared_gene_stats'].keys()),
        'Value': list(stats['overlap_statistics'].values()) + list(stats['shared_gene_stats'].values())
    })
    overlap.to_csv(f"{output_prefix}_overlap_stats.tsv", sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Compare GTF files and generate statistics')
    parser.add_argument('--reference', required=True, help='Reference GTF file')
    parser.add_argument('--query', required=True, help='Query GTF file to compare')
    parser.add_argument('--output-prefix', required=True, help='Prefix for output files')
    
    args = parser.parse_args()
    
    try:
        stats = compare_gtfs(args.reference, args.query)
        write_stats(stats, args.output_prefix)
        logger.info(f"Statistics written to files with prefix: {args.output_prefix}")
    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}")
        raise

if __name__ == '__main__':
    main()