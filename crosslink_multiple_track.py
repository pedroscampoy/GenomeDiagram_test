#!/home/pjsola/env/bin/python

import os
import csv

from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink


diagram_name = 'TEST_CL'
gd_diagram = GenomeDiagram.Diagram(diagram_name)
dict_records = {'NC_016838.1':122799, 'NC_016839.1':105974, 'NC_016846.1':111195}

#NC_016838.1 vs NC_016839.1 made up reltions
A_vs_B = [
    (99, "mcpQ", "tetR"),
    (33, "04374", "xerD")
]

A_vs_B = [
    (99, "tetA", "pld"),
    (33, "rhsC", "traC")
]
i = 0
for record,record_length in dict_records.items():
    # Allocate tracks 5 (top), 3, 1 (bottom) for A, B, C
    # (empty tracks 2 and 4 add useful white space to emphasise the cross links
    # and also serve to make the tracks vertically more compressed)
    

    gd_track_for_features = gd_diagram.new_track(5 - 2 * i,name=record, greytrack=True, height=0.5,
                                                 start=0, end=record_length)

    gd_set_features = gd_track_for_features.new_set()
    #gd_track_for_features = gd_diagram.new_track(1, name=record, greytrack=True, start=0, end=record_length)
    #gd_set_features = gd_track_for_features.new_set()

    i += 1


    with open('KPN.gff.forward.coordinates', 'r') as bed_forward_file:
        bed_readed = csv.reader(bed_forward_file, delimiter="\t")

        #record = None
        max_position = 0

        for row in bed_readed:
            #NC_016839.1	122155	122652	04281
            #record  loc1    loc2    name
            if row[0] == record:
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
            #NC_016839.1	122155	122652	04281
            #record  loc1    loc2    name
            if row[0] == record:
                #if record is None:
                #    record = 'row[0]'
            #Add three features to show the strand options,
                if int(row[2]) > max_position:
                    max_position = int(row[2]) 

                feature = SeqFeature(FeatureLocation(int(row[1]), int(row[2])), strand=-1)
                gd_set_features.add_feature(feature, name=row[3], label=True, color="green",
                            label_size=6, label_angle=90, label_position="middle", sigil="ARROW", arrowshaft_height=1.0)


#for score, x, y in A_vs_B:
#color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick, 0, 100, 50)
#    border = colors.lightgrey 	
#link_xy = CrossLink((5, 321, 1418), (3, 5544, 6503), color=colors.lightgrey)
#gd_diagram.cross_track_links.append(link_xy)

print('TEST')
print(gd_set_features.id)



# track_X = gd_diagram.tracks[5]
# track_Y = gd_diagram.tracks[3]
#     for score, id_X, id_Y in A_vs_B:
#         feature_X = get_feature(rec_X.features, id_X)
#         feature_Y = get_feature(rec_Y.features, id_Y)
#         color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick, 0, 100, score)
#         link_xy = CrossLink((track_X, feature_X.location.start, feature_X.location.end),
#                             (track_Y, feature_Y.location.start, feature_Y.location.end),
#                             color, colors.lightgrey)
#         gd_diagram.cross_track_links.append(link_xy)

#gd_diagram.draw(format='linear', pagesize=(15*cm,4*cm), fragments=1,start=0, end=max_position)
#gd_diagram.draw(format='linear', pagesize='A4', fragments=1,start=0, end=max_position)

#gd_diagram.write( diagram_name + ".pdf", "pdf")