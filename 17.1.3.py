#!/home/pjsola/env/bin/python

import os
from Bio import Entrez
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO


#os.remove("NC_005816.gb")
if os.path.exists("NC_005816.gb"):
    print('NC_005816.gb already exist')
else:
    print('NC_005816.gb DO NOT exist')
    handle = Entrez.efetch(db="nucleotide", id="NC_005816", rettype="gb", retmode="text")
    with open("NC_005816.gb","w") as file_genebank:
        #file_genebank = open("NC_005816.gb","w")
        file_genebank.write(handle.read())
        #file_genebank.close
 
record = SeqIO.read("NC_005816.gb", "genbank")
#print(record.__dict__)

#after loading in our sequence we next create an empty diagram, then add an (empty) track, and to that add an (empty) feature set
gd_diagram = GenomeDiagram.Diagram(record.description)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type != "gene":
        #Exclude this feature
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
        #print (feature.__dict__)
        #print (feature.location)
        #print (feature.qualifiers['locus_tag'])
        #print (feature.qualifiers['old_locus_tag'])
    else:
        color = colors.lightblue

    gd_feature_set.add_feature(feature, color=color, label=True)

gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4', fragments=1, start=0, end=len(record))
#gd_diagram.write("plasmid_linear.pdf", "PDF")
#gd_diagram.write("plasmid_linear.eps", "EPS")
#gd_diagram.write("plasmid_linear.svg", "SVG")
gd_diagram.write("plasmid_linear.png", "PNG")
#gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm), start=0, end=len(record), circle_core=0.7)
