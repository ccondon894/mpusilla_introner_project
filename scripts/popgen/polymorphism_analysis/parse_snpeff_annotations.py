#!/usr/bin/env python3

import pandas as pd
import re
import sys

def parse_snpeff_annotations():
    """
    Parse SnpEff annotations from the INFO field to classify mutation impacts
    """
    print("=== SnpEff Annotation Parsing ===")
    print("Step 4: Parsing mutation impact classifications")
    
    # Read the extracted variants
    variants_df = pd.read_csv('/scratch1/chris/introner_vis/introner_variants.tsv', sep='\t')
    print(f"Processing {len(variants_df)} variants")
    
    # Parse ANN field from INFO column
    annotation_data = []
    
    for idx, row in variants_df.iterrows():
        if idx % 1000 == 0:
            print(f"Processed {idx}/{len(variants_df)} variants")
            
        info_field = str(row['info'])
        ortholog_id = row['ortholog_id']
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        
        # Extract ANN field (SnpEff annotations)
        ann_match = re.search(r'ANN=([^;]+)', info_field)
        if ann_match:
            ann_field = ann_match.group(1)
            
            # Multiple annotations are separated by commas
            annotations = ann_field.split(',')
            
            for ann in annotations:
                # Parse individual annotation
                # Format: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos|CDS.pos|Protein.pos|Distance|Warnings
                fields = ann.split('|')
                
                if len(fields) >= 15:
                    allele = fields[0]
                    annotation = fields[1]
                    impact = fields[2]
                    gene_name = fields[3]
                    gene_id = fields[4]
                    feature_type = fields[5]
                    feature_id = fields[6]
                    transcript_biotype = fields[7]
                    rank = fields[8]
                    hgvs_c = fields[9]
                    hgvs_p = fields[10]
                    cdna_pos = fields[11]
                    cds_pos = fields[12]
                    protein_pos = fields[13]
                    distance = fields[14]
                    
                    # Classify mutation type
                    mutation_type = classify_mutation_type(annotation, impact)
                    
                    annotation_data.append({
                        'ortholog_id': ortholog_id,
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'allele': allele,
                        'annotation': annotation,
                        'impact': impact,
                        'mutation_type': mutation_type,
                        'gene_name': gene_name,
                        'gene_id': gene_id,
                        'feature_type': feature_type,
                        'transcript_biotype': transcript_biotype,
                        'hgvs_c': hgvs_c,
                        'hgvs_p': hgvs_p
                    })
        else:
            # No ANN field found
            annotation_data.append({
                'ortholog_id': ortholog_id,
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'allele': alt,
                'annotation': 'no_annotation',
                'impact': 'unknown',
                'mutation_type': 'unknown',
                'gene_name': '',
                'gene_id': '',
                'feature_type': '',
                'transcript_biotype': '',
                'hgvs_c': '',
                'hgvs_p': ''
            })
    
    # Convert to DataFrame
    annotations_df = pd.DataFrame(annotation_data)
    
    print(f"\nAnnotation results:")
    print(f"Total annotations: {len(annotations_df)}")
    print(f"Unique variants: {annotations_df[['chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]}")
    
    # Summarize mutation types
    print("\nMutation type distribution:")
    print(annotations_df['mutation_type'].value_counts())
    
    print("\nAnnotation type distribution:")
    print(annotations_df['annotation'].value_counts().head(10))
    
    print("\nImpact distribution:")
    print(annotations_df['impact'].value_counts())
    
    # Save annotations
    output_file = '/scratch1/chris/introner_vis/variant_annotations.tsv'
    annotations_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nAnnotations saved to: {output_file}")
    
    return annotations_df

def classify_mutation_type(annotation, impact):
    """
    Classify mutation into synonymous, missense, nonsense, or other
    """
    annotation = annotation.lower()
    impact = impact.lower()
    
    # Synonymous variants
    if 'synonymous' in annotation:
        return 'synonymous'
    
    # Missense variants
    if 'missense' in annotation:
        return 'missense'
    
    # Nonsense/stop variants
    if any(x in annotation for x in ['stop_gained', 'nonsense', 'stop_lost']):
        return 'nonsense'
    
    # Start lost
    if 'start_lost' in annotation:
        return 'start_lost'
    
    # Frameshift
    if 'frameshift' in annotation:
        return 'frameshift'
    
    # Splice site variants
    if any(x in annotation for x in ['splice_donor', 'splice_acceptor', 'splice_region']):
        return 'splice'
    
    # Intron variants
    if 'intron' in annotation:
        return 'intron'
    
    # UTR variants
    if any(x in annotation for x in ['3_prime_utr', '5_prime_utr', 'utr']):
        return 'utr'
    
    # Upstream/downstream
    if any(x in annotation for x in ['upstream', 'downstream']):
        return 'regulatory'
    
    # Intergenic
    if 'intergenic' in annotation:
        return 'intergenic'
    
    # Default to impact level if specific annotation not recognized
    if impact in ['high', 'moderate', 'low', 'modifier']:
        return f'impact_{impact}'
    
    return 'other'

if __name__ == "__main__":
    result = parse_snpeff_annotations()
    if result is not None:
        print(f"\n✓ Successfully parsed {len(result)} annotations")
    else:
        print("\n✗ Failed to parse annotations")
        sys.exit(1)