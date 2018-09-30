#!/home/pjsola/env/bin/python

import os
import csv

from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

gdd = GenomeDiagram.Diagram('NC_016838.1')
gd_track_for_features = gdd.new_track(1, greytrack=True)
gd_set_features = gd_track_for_features.new_set()


with open('KPN.gff.forward.coordinates', 'r') as bed_forward_file:
    bed_readed = csv.reader(bed_forward_file, delimiter="\t")

    #record = None
    max_position = 0

    for row in bed_readed:
        #NC_016838.1	122155	122652	04281
        #record  loc1    loc2    name
        if row[0] == 'NC_016838.1':
            #if record is None:
            #    record = 'row[0]'
        #Add three features to show the strand options,
            if int(row[2]) > max_position:
                max_position = int(row[2]) 

            feature = SeqFeature(FeatureLocation(int(row[1]), int(row[2])), strand=+1)
            gd_set_features.add_feature(feature, name=row[3], label=True, color="blue",
                           label_size=6, label_angle=90, label_position="middle", sigil="ARROW", arrowshaft_height=1.0)


with open('KPN.gff.reverse.coordinates', 'r') as bed_reverse_file:
    bed_readed_rev = csv.reader(bed_reverse_file, delimiter="\t")

    #record = None

    for row in bed_readed_rev:
        #NC_016838.1	122155	122652	04281
        #record  loc1    loc2    name
        if row[0] == 'NC_016838.1':
            #if record is None:
            #    record = 'row[0]'
        #Add three features to show the strand options,
            if int(row[2]) > max_position:
                max_position = int(row[2]) 

            feature = SeqFeature(FeatureLocation(int(row[1]), int(row[2])), strand=-1)
            gd_set_features.add_feature(feature, name=row[3], label=True, color="green",
                           label_size=6, label_angle=90, label_position="middle", sigil="ARROW", arrowshaft_height=1.0)


#gdd.draw(format='linear', pagesize=(15*cm,4*cm), fragments=1,start=0, end=max_position)
gdd.draw(format='linear', pagesize='A4', fragments=4,start=0, end=max_position)

gdd.write("NC_016838.1_all_arrrow_grey.pdf", "pdf")