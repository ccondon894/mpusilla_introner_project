import sys
import requests
import json
import pandas as pd

tsv = sys.argv[1]
df = pd.read_csv(tsv, sep="\t")

outfile = sys.argv[2]

ko_dict = {}
for index, row in df.iterrows():
    if pd.notna(row['KO']):
        ko_ids = row['KO'].split(',')
        for id in ko_ids:
            ko_dict[id] = set()


def get_go_terms_From_ko(ko_id):
    url = f"https://rest.kegg.jp/get/{ko_id}"
    response = requests.get(url)
    go_terms = set()
    
    if response.status_code == 200:
        lines = response.text.strip().split('\n')
        in_dblinks_section = False
        
        for line in lines:
            if line.startswith("DBLINKS"):
                in_dblinks_section = True
                # Process this line for GO terms
                process_line_for_go_terms(line, go_terms)
            elif in_dblinks_section and line.startswith("            "):
                # Process continuation lines that are part of DBLINKS section
                process_line_for_go_terms(line, go_terms)
            elif in_dblinks_section and not line.startswith("            "):
                # We've reached the end of the DBLINKS section
                in_dblinks_section = False
    
    return go_terms

def process_line_for_go_terms(line, go_terms):
    parts = line.strip().split()
    for i, part in enumerate(parts):
        if part == "GO:":
            go_terms.update([f"GO:{go}" for go in parts[i+1:]])
        elif part == "GO" and i+1 < len(parts) and parts[i+1] == ":":
            go_terms.update([f"GO:{go}" for go in parts[i+2:]])
        elif part.startswith("GO:"):
            go_terms.add(part)


for id in ko_dict.keys():
    print(f"Querying id {id}")
    go_terms = get_go_terms_From_ko(id)
    ko_dict[id].update(go_terms)

ko_dict_for_json = {}

for key, value in ko_dict.items():
    ko_dict_for_json[key] = list(value)

with open(outfile, 'w') as f:
    json.dump(ko_dict_for_json, f, indent=2)
