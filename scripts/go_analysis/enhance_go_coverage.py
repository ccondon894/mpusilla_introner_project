import pandas as pd
import requests
import json
import sys
from collections import defaultdict

def flatten_mixed_list(mixed_list):
    """Recursively flatten a mixed list of strings and lists."""
    result = []
    for item in mixed_list:
        if isinstance(item, list):
            result.extend(flatten_mixed_list(item))  # Recursive call for nested lists
        elif isinstance(item, tuple):
            result.extend(flatten_mixed_list(list(item)))  # Handle tuples as well
        else:
            result.append(item)  # Add the string directly
    return result

def process_file(filename):
    """Safely load a JSON file."""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return {}

def main():
    # Process command line arguments
    if len(sys.argv) < 7:
        print("Usage: script.py tsv pfam_json ko_json tair_json panther_json outfile")
        sys.exit(1)
        
    tsv = sys.argv[1]
    pfam_json = sys.argv[2]
    ko_json = sys.argv[3]
    tair_json = sys.argv[4]
    panther_json = sys.argv[5]
    outfile = sys.argv[6]
    
    # Load the TSV file
    print("Loading annotation file...")
    df = pd.read_csv(tsv, sep="\t")
    
    # Initialize gene_to_go with defaultdict for simpler code
    gene_to_go = defaultdict(set)
    
    # Parse existing GO terms
    print("Processing existing GO terms...")
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['GO']):
            gene_to_go[gene_id].update(row['GO'].split(','))
    
    # Process PFAM terms
    print("Processing Pfam terms...")
    pfam_data = process_file(pfam_json)
    
    # Convert to sets for efficiency
    pfam_dict = {}
    for key, lists in pfam_data.items():
        pfam_dict[key] = set(item for sublist in lists for item in sublist)  # Flatten and convert to set
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['Pfam']):
            for pfam_id in row['Pfam'].split(','):
                if pfam_id in pfam_dict:
                    gene_to_go[gene_id].update(pfam_dict[pfam_id])
    
    # Process KO terms
    print("Processing KO terms...")
    ko_dict = process_file(ko_json)
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['KO']):
            ko_id = row['KO']
            if ko_id in ko_dict:
                gene_to_go[gene_id].update(ko_dict[ko_id])
    
    # Process TAIR terms
    print("Processing TAIR terms...")
    tair_dict = process_file(tair_json)
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['Best-hit-arabi-name']):
            tair_id = row['Best-hit-arabi-name']
            if '.' in tair_id:
                tair_id = tair_id.split('.')[0]
            if tair_id in tair_dict:
                gene_to_go[gene_id].update(tair_dict[tair_id])
    
    # Process PANTHER json
    print("Processing PANTHER terms...")
    panther_dict = process_file(panther_json)
    
    for index, row in df.iterrows():
        gene_id = row['locusName']
        if pd.notna(row['Panther']):
            for panther_id in row['Panther'].split(','):
                if panther_id in panther_dict:
                    gene_to_go[gene_id].update(panther_dict[panther_id])
    
    # Prepare final output
    print("Preparing output...")
    gene_to_go_for_json = {}
    
    for gene, terms in gene_to_go.items():
        # Convert to list, flatten, and ensure uniqueness
        flattened_terms = flatten_mixed_list(list(terms))
        gene_to_go_for_json[gene] = list(set(flattened_terms))
    
    # Write the output
    print(f"Writing results to {outfile}...")
    with open(outfile, 'w') as f:
        json.dump(gene_to_go_for_json, f, indent=2)
    
    print("Done!")

if __name__ == "__main__":
    main()