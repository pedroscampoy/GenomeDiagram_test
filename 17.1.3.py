#!/usr/bin/python3

from Bio import Entrez

handle = Entrez.efetch(db="nucleotide", id="NC_005816", rettype="gb", retmode="text")


with open("NC_005816.gb","w") as file_genebank:

#file_genebank = open("NC_005816.gb","w")
    file_genebank.write(handle.read())
#file_genebank.close
