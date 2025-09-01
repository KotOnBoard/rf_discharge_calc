import json
import pandas as pd
cname = 'Ar_IzCross'
with open(f"{cname}.json5", 'r') as file:
    cross = file.readlines()
    for i in range(0,len(cross)):
        cross[i] = pd.to_numeric(cross[i].split())
