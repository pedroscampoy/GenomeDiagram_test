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
dict_records = {'NC_016838':122799, 'NC_016839':105974, 'NC_016846':111195}

#NC_016838.1 vs NC_016839.1 made up reltions

features_x_link = []
features_y_link = []

B_vs_C = [
    (99, "tetA", "pld"),
    (33, "rhsC", "traC")
]
i = 0

for record,record_length in dict_records.items():
    # Allocate tracks 5 (top), 3, 1 (bottom) for A, B, C
    # (empty tracks 2 and 4 add useful white space to emphasise the cross links
    # and also serve to make the tracks vertically more compressed)
    

    gd_track_for_features = gd_diagram.new_track(5 - 2 * i,name=record, greytrack=True, greytrack_labels=5, scale_ticks=0,height=0.5,
                                                 start=0, end=record_length)
    name_for_featureset = 'gd_set_features_'+record
    print('--------------------------------------------------------------->'+name_for_featureset)

    name_for_featureset = gd_track_for_features.new_set(name=record)
    #gd_track_for_features = gd_diagram.new_track(1, name=record, greytrack=True, start=0, end=record_length)
    #gd_set_features = gd_track_for_features.new_set()

    
    count_feature = 0
    num = 0

    with open('KPN.gff.forward.coordinates', 'r') as bed_forward_file:
        bed_readed = csv.reader(bed_forward_file, delimiter="\t")

        #record = None

        for row in bed_readed:
            #NC_016839.1	122155	122652	04281
            #record  loc1    loc2    name
            if row[0] == record:
                
                #if record is None:
                #    record = 'row[0]'
                #Add three features to show the strand options,
                #print('NUM IS ' + str(num))

                feature = SeqFeature(FeatureLocation(int(row[1]), int(row[2])), strand=+1)
                name_for_featureset.add_feature(feature, name=row[3], label=True, color="blue",
                            label_size=6, label_angle=90, label_position="middle", sigil="ARROW", arrowshaft_height=1.0)

                #print(name_for_featureset[int(num)].name)
                num += 1


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
                
                #print('NUM IS ' + str(num))


                feature = SeqFeature(FeatureLocation(int(row[1]), int(row[2])), strand=-1)
                name_for_featureset.add_feature(feature, name=row[3], label=True, color="green",
                            label_size=6, label_angle=90, label_position="middle", sigil="ARROW", arrowshaft_height=1.0)

                #print(name_for_featureset[int(num)].name)
                num += 1
    
    A_vs_B = [
    (50, "mcpQ", "tetR"),
    (33, "uvrB", "rhsC"),
    (10, "04321", "05658"),
    (100, "ligA", "soj"),
    (100, "ligA", "korB"),
    (100, "smc", "xerD")
]
    

    for number, cross_link_relation in enumerate(A_vs_B):
        score = cross_link_relation[0]
        cross_link_feature_A = cross_link_relation[1]
        cross_link_feature_B = cross_link_relation[2]
        print ("NUMBER IS:" + str(number))

        for feature_number in range(1,len(name_for_featureset)):

            if name_for_featureset[feature_number].name == cross_link_feature_A:
                
                features_x_link.insert(number, name_for_featureset[feature_number])

                print ("Test FEATUREA 1: ")
                print (number)
                print (name_for_featureset[feature_number].name)
                print(features_x_link)

                track_x_name = name_for_featureset.name
                track_x = name_for_featureset
                #print ("Test FEATUREA: " + features_x_link[number])
                #print (number)
                


            if name_for_featureset[feature_number].name == cross_link_feature_B and name_for_featureset.name != str(track_x_name):
                features_y_link.insert(number, name_for_featureset[feature_number])
                #feature_y_n = name_for_featureset[feature_number]
                track_y_name = name_for_featureset.name
                track_y = name_for_featureset

                print ("Test FEATUREA 2: ")
                print (number)
                print (name_for_featureset[feature_number].name)
                print (features_y_link)

               
    

    for feature_x, feature_y in zip(features_x_link, features_y_link):
        #print("SCORE: " + str(count_feature))
        score = A_vs_B[count_feature][0]
        color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick,
                                             0, 100, score)
        border = colors.lightgrey
        gd_diagram.cross_track_links.append(CrossLink(feature_x, feature_y, color, border))
        count_feature += 1

    i += 1
    

gd_diagram.draw(format='linear', pagesize='A4', fragments=1,start=0)

gd_diagram.write( diagram_name + ".pdf", "pdf")
gd_diagram.write( diagram_name + ".svg", "SVG")