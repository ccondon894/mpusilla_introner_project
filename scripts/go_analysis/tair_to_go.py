import sys
import pandas as pd
import json

tsv = sys.argv[1]
tair_assoc = sys.argv[2]
outfile = sys.argv[3]

df = pd.read_csv(tsv, sep="\t")
tair_dict = {}
for index, row in df.iterrows():
    if pd.notna(row['Best-hit-arabi-name']):
        id = row['Best-hit-arabi-name']
        if '.' in id:
            id = id.split('.')[0]
        
        tair_dict[id] = set()


with open(tair_assoc, 'r') as f:
    for line in f:
        if line.startswith('!'):
            continue

        info = line.strip().split("\t")
        id = info[1]
        go_term = info[4]

        if id in tair_dict:
            tair_dict[id].add(go_term)

tair_dict_for_json = {}
for key, value in tair_dict.items():
    tair_dict_for_json[key] = list(value)

with open(outfile, 'w') as f:
    json.dump(tair_dict_for_json, f, indent=2)