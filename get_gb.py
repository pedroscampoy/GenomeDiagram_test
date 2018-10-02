#!/home/pjsola/env/bin/python

import os
from Bio import Entrez

handle = Entrez.efetch(db="nucleotide", id="NC_003212", rettype="gb", retmode="text")
with open("NC_003212.gbk","w") as file_genebank:
    file_genebank.write(handle.read())