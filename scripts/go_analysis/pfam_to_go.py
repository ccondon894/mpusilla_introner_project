import sys
import re
import json
import requests
import pandas as pd
from collections import defaultdict

def read_interpro2go(interpro2go):

    interpro2go_dict = defaultdict(list)
    go_pattern = r'^GO:\d{7}$'
    with open(interpro2go, 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue

            info = line.strip().split(' ')
            
            if not info[0].startswith('InterPro'):
                continue
            ipr_id = info[0].split(":")[1]
            for item in info[1:]:
                if bool(re.match(go_pattern, item)):
                    interpro2go_dict[ipr_id].append(item)

    return interpro2go_dict


def get_interpro_from_pfam(pfam_id):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}"
    response = requests.get(url, headers={"Accept": "application/json"})
    
    if response.status_code == 200:
        data = response.json()
        
        # Check if go_terms exists in metadata and is not null
        metadata = data.get('metadata', {})
        ipr_id = metadata.get("integrated")
        return ipr_id
    else:
        return None

interpro2go = sys.argv[1]
tsv = sys.argv[2]
outfile = sys.argv[3]

my_dict = read_interpro2go(interpro2go)
my_dict = dict(my_dict)


df = pd.read_csv(tsv, sep="\t")
pfam_to_go_dict = defaultdict(list)
for index, row in df.iterrows():
    if pd.notna(row['Pfam']):
        pfam_ids = row['Pfam'].split(',')
        for id in pfam_ids:
            ipr_id = get_interpro_from_pfam(id)
            if ipr_id and ipr_id in my_dict:
                go_terms = my_dict[ipr_id]
                pfam_to_go_dict[id].append(go_terms)
            else:
                print("uh oh")


pfam_to_go_dict = dict(pfam_to_go_dict)
with open(outfile, 'w') as f:
    json.dump(pfam_to_go_dict, f, indent=2)