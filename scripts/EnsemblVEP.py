import requests, sys
import pandas as pd

server = "https://grch37.rest.ensembl.org"
ext = "/vep/human/hgvs"
headers={ "Content-Type" : "application/json", "Accept" : "application/json" }
r = requests.post(server+ext, headers=headers, data='{ "hgvs_notations" : ["9:g.22125503G>C"], "CADD": "1", "GeneSplicer": "1", "Phenotypes": "1", "canonical": "1", "dbscSNV": "1", "domains": "1", "hgvs": "1", "merged": "1", "miRNA": "1", "numbers": "1", "protein": "1", "refseq": "1", "uniprot": "1", "variant_class": "1", "xref_refseq": "1" }')
 
if not r.ok:
	r.raise_for_status()
	sys.exit()
 
decoded = r.json()
print(decoded[0])

