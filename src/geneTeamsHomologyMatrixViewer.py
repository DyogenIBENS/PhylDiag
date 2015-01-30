#!/usr/bin/python
# -*- coding: utf-8 -*-

# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import os
import sys
import itertools

#import utils.myGenomes as myGenomes
import utils.myLightGenomes as myLightGenomes
import utils.myTools as myTools
# PhylDiag core algorithm
import utils.myDiags as myDiags
import utils.myGeneTeams as myGeneTeams

import drawHomologyMatrixWithSBs

__doc__= """
        Show the homology matrix with coloured gene teams.
        - On the x-axis are the genes of the 1st genome in the desired window
        - On the y-axis are the genes of the 2nd genome in the desired window
        - Each coloured rectangle in the matrix represents a homology
        - '+' indicates homology that have horizontal gene and vertical gene in the same direction on both chromosomes
          '-' indicates homology that have horizontal gene and vertical gene in opposite directions on each chromosomes
        """


def genesComputeGeneTeamsIndices(listOfGeneTeams):
    ###
    # Build geneTeamsIndices = [..., [...,(i16,j16),...], ...] list of gene
    # teams with
    # geneteams = list of all the points of the diagonal
    ###
    gtIndices = []
    for gt in listOfGeneTeams:
        (la, (c1, l1), (c2, l2)) = gt
        gtIndices.append([])
        for (idxAG, _) in enumerate(la):
            li1s = [i1 for i1s in l1[idxAG] for i1 in i1s]
            li2s = [i2 for i2s in l2[idxAG] for i2 in i2s]
            gtIndices[-1].extend(itertools.product(li1s, li2s))
    return gtIndices

arguments = myTools.checkArgs(
    [("genome1", file), ("genome2", file),
     ("families", file),
     ("chr1:deb1-fin1", str),
     ("chr2:deb2-fin2", str)],
    [("gapMax", str, 'None'),
     ("tandemGapMax", int, 0),
     ("filterType", str, 'InFamilies'),
     ("minChromLength", int, 1),
     ("out:GeneTeams", str, "./res/geneTeamsDrawer.txt"),
     #TODO ("mode:chromosomesRewrittenInTbs", bool, False),
     #TODO ('convertGenicToTbCoordinates', bool, False),
     ("out:ImageName", str, "./res/homologyMatrix.svg"),
     ('verbose', bool, True)],
    __doc__)

#by convention:
if arguments['gapMax'] == 'None':
    arguments['gapMax'] = None
else:
    try:
        arguments['gapMax'] = int(arguments['gapMax'])
    except:
        raise ValueError('gapMax must be an int or None')

#assert (arguments["convertGenicToTbCoordinates"] and
#        arguments["mode:chromosomesRewrittenInTbs"])\
#        or not arguments["convertGenicToTbCoordinates"]

# Load genomes
genome1 = myLightGenomes.LightGenome(arguments["genome1"])
genome2 = myLightGenomes.LightGenome(arguments["genome2"])
# Change genome format
genome1Name = genome1.name
genome2Name = genome2.name
genome1_ = {}
families = myLightGenomes.Families(arguments["families"])
modesFilter = list(myDiags.FilterType._keys)
filterType = myDiags.FilterType[modesFilter.index(arguments["filterType"])]
#thresholdChr = 50
#chromosomes are shown as a list of genes
#print >> sys.stderr, "List of (chromosomes, length in genes) of Genome 1, for chr of size > %s " % thresholdChr
#for (chr1,len1) in [(key1, len(chr1)) for (key1,chr1) in genome1.items() if len(chr1) > thresholdChr]:
#       print >> sys.stderr, "chr %s has %s genes" % (chr1,len1)
#print >> sys.stderr, "List of (chromosomes, length in genes) of Genome 2, for chr of size > %s " % thresholdChr
#for (chr2,len2) in [(key2, len(chr2)) for (key2,chr2) in genome2.items() if len(chr2) > thresholdChr]:
#       print >> sys.stderr, "chr %s has %s genes" % (chr2,len2)
(chr1,range1) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr1:deb1-fin1"], genome1)
(chr2,range2) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr2:deb2-fin2"], genome2)
chrom1 = myLightGenomes.LightGenome()
chrom2 = myLightGenomes.LightGenome()
chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]

###
# Build Genes Strands
###
genesStrandsC1 = [s for (_,s) in chrom1[chr1]]
genesStrandsC2 = [s for (_,s) in chrom2[chr2]]

((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
 genesHomologiesHpSign,
 (genesNoHomologiesInWindowC1,genesNoHomologiesInWindowC2),
 genesHomologyGroupsInWindow) =\
    drawHomologyMatrixWithSBs.genesComputeHomologyInformations(chr1, chr2, chrom1, chrom2,
                                          families,
                                          filterType,
                                          arguments['minChromLength'],
                                          arguments['tandemGapMax'])

# Search gene teams
#FIXME : calculate gene teams before, on the whole chromosome, not the ROI specified by the user ranges
listOfGeneTeams =\
    list(myGeneTeams.extractGtsInPairCompGenomes(chrom1,
                                                 chrom2,
                                                 families,
                                                 tandemGapMax=arguments['tandemGapMax'],
                                                 gapMax=arguments['gapMax'],
                                                 filterType=filterType,
                                                 minChromLength=arguments["minChromLength"],
                                                 verbose=arguments['verbose']))
genesGeneTeamsIndices = genesComputeGeneTeamsIndices(listOfGeneTeams)

strArray = drawHomologyMatrixWithSBs.drawHomologyMatrix(
    (range1, range2), (genesStrandsC1, genesStrandsC2),
    (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
    (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
    genesHomologiesHpSign, genesHomologyGroupsInWindow,
    genesGeneTeamsIndices,
    outputFileName=arguments["out:ImageName"], maxWidth=100, maxHeight=100 )

#copy the css style sheet
dirNameImage = os.path.dirname(arguments["out:ImageName"])
dirNameImage = dirNameImage if dirNameImage != "" else "."
print >> sys.stderr, "cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage)
os.system("cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage))

# write a simple file with all diagonals into output file
f=open(arguments['out:GeneTeams'],'w')
print >> f, "Mode : %s" %  'Genic scale'
print >> f, "chromosome %s de %s\t%s\t%s\tchromosome %s de %s\t%s\t%s\t%s" % (chr1, genome1Name, 'beginC1', 'endC1', chr2, genome2Name, 'beginC2', 'endC2','length in ancestral genes')
print >> f, "c1\tbeg1\tend1\tc2\tbeg2\tend2\thps\tpVal"
listOfGeneTeams = [l for l in listOfGeneTeams]
listOfGeneTeams = sorted(list(listOfGeneTeams), key=lambda x:len(x[2]), reverse=True)
for gt in listOfGeneTeams:
    (la, (c1,l1), (c2,l2)) = gt
    min_l2 = min(l2[0],l2[-1])
    max_l2 = max(l2[0],l2[-1])
    print >> f, "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (c1, l1[0], l1[-1], c2, min_l2, max_l2, len(la))
f.close()

# Add lengends and title to the ouput matrix
height=100
width=100
var = ['<?xml version="1.0" encoding="utf-8" standalone="no"?>\n',
                '<?xml-stylesheet type="text/css" href="styleForHomologyMatrixWithSBs.css" ?>\n', #Warning : request the css file
                '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n',
                "<svg height=\"100%%\" version=\"1.1\" viewBox=\"0 0 %s %s\" width=\"100%%\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" % (width, height),
                 '<defs>\n',
                  '<style type="text/css">\n',
                        '*{stroke-linecap:square;stroke-linejoin:round;}\n',
                  '</style>\n',
                 '</defs>\n'
                 '<g style="fill-opacity:1.0; stroke:black;\n',
                 'stroke-width:1;">\n']
#Title
title =\
    "%s, tandemGapMax=%s tbs, gapMax=%s tbs, %s sbs" %\
    ('MH',
     arguments['tandemGapMax'],
     arguments['gapMax'],
     len(listOfGeneTeams))

var += ['<svg x="5" y="0" viewBox="5 0 95 5" width="95" height="5" xmlns="http://www.w3.org/2000/svg" xmlns:xhtml="http://www.w3.org/1999/xhtml" xmlns:xlink="http://www.w3.org/1999/xlink">\n',
                '<foreignObject x="0" y="0" width="95" height="5">\n',
                '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
                        '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
                        '<xhtml:div style="color:black; word-wrap:break-word; font-size:2px; font-family:Arial" >' + title + '\n',
                                        '</xhtml:div>\n',
                                '</xhtml:div>\n',
                        '</xhtml:div>\n',
                '</foreignObject>\n',
        '</svg>\n']

#Add legends (genomes names and ranges)
var += ['<svg x="0" y="0" viewBox="0 0 5 95" width="5" height="95" xmlns="http://www.w3.org/2000/svg" xmlns:xhtml="http://www.w3.org/1999/xhtml" xmlns:xlink="http://www.w3.org/1999/xlink">\n',
                '<foreignObject x="0" y="0" width="95" height="5" transform="translate(5,0) rotate(90) translate(0,0)">\n',
                '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
                        '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
                        '<xhtml:div style="color:black; word-wrap:break-word; font-size:2px; font-family:Arial" >' + str(arguments["genome2"]) + " <br />  chr" + str(chr2) + ":" + str(range2[0]+1) + "-" + str(range2[1]) + '\n',
                                        '</xhtml:div>\n',
                                '</xhtml:div>\n',
                        '</xhtml:div>\n',
                '</foreignObject>\n',
        '</svg>\n',
        '<svg x="5" y="95" viewBox="0 0 95 5" width="95" height="5" xmlns="http://www.w3.org/2000/svg" xmlns:xhtml="http://www.w3.org/1999/xhtml" xmlns:xlink="http://www.w3.org/1999/xlink">\n',
                '<foreignObject x="0" y="0" width="95" height="5">\n',
                '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
                        '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
                        '<xhtml:div style="color:black; word-wrap:break-word; font-size:2px; font-family:Arial" >' + str(arguments["genome1"]) + " <br /> chr" + str(chr1) + ":" + str(range1[0]+1) + "-" + str(range1[1]) + '\n',
                                        '</xhtml:div>\n',
                                '</xhtml:div>\n',
                        '</xhtml:div>\n',
                '</foreignObject>\n',
        '</svg>\n']

var+=['<svg preserveAspectRatio="xMidYMid meet" x="5" y="5" viewBox="0 0 100 100" width="90" height="90" >\n'] # little transformation : viewBox = "the part of the forecoming images that we want to see", width and height = the width and height of the image that will be printed on the screen. This instructions takes a viewBox of the forecoming images and display it in a image of the specified width and height
for line in strArray:
    if line.find("<?xml")>=0 or line.find("<!DOCTYPE")>=0: #or line.find("<svg")>=0 or line.find("</svg")>=0:
        continue
    else:
        var += line

var += ["</svg>\n"]
var += ["</g>\n", "</svg>\n"]

file = open(arguments["out:ImageName"],'w')
file.writelines(var)
file.close()
#os.system("%s %s" % ('firefox',arguments["out:ImageName"]))
