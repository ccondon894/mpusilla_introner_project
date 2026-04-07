import sys
import pandas as pd
import json

panther_db = sys.argv[1]
outfile = sys.argv[2]

# Process PANTHER terms
panther_dict = {}
with open(panther_db, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) >= 2:
            panther_id = fields[0]
            go_terms = []
            for field in fields[3:]:
                if 'GO:' in field:
                    term_chunks = field.split(';')
        
                    for chunk in term_chunks:
                        # Extract the GO ID using the pattern #GO:
                        if '#GO:' in chunk:
                            go_id = chunk.split('#')[1]
                            go_terms.append(go_id)
        
        panther_dict[panther_id] = set(go_terms)


panther_dict_for_json = {}
for key, values in panther_dict.items():
    panther_dict_for_json[key] = list(values)

with open(outfile, 'w') as f:
    json.dump(panther_dict_for_json, f, indent=2)


