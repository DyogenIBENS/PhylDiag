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

import utils.myLightGenomes as myLightGenomes
import utils.myTools as myTools
# PhylDiag core algorithm
import utils.myDiags as myDiags
import utils.myMapping as myMapping

import drawHomologyMatrixWithSBs


__doc__ = """
        Show the homology matrix with coloured synteny blocks (also called diagonals).
        - On the x-axis are the genes of the 1st genome in the desired window
        - On the y-axis are the genes of the 2nd genome in the desired window
        - Each coloured rectangle in the matrix represents a homology
        - '+' indicates homology that have horizontal gene and vertical gene in the same direction on both chromosomes
          '-' indicates homology that have horizontal gene and vertical gene in opposite directions on each chromosomes
        """

if __name__ == '__main__':

    arguments = myTools.checkArgs(
        [("genome1", file),
         ("genome2", file),
         ("families", file),
         ("chr1:deb1-fin1", str),
         ("chr2:deb2-fin2", str)],
        [("filterType", str, 'InFamilies'),
         ("tandemGapMax", int, 0),
         ("gapMax", str, 'None'),
         ('identifyBreakpointsWithinGaps', bool, False),
         ("nonOverlappingSbs", bool, False),
         ("overlapMax", int, 0),
         ("consistentSwDType", bool, True),
         ("considerMonogenicSb", bool, False),
         ("minChromLength", int, 1),
         ("pThreshold", float, 0.001),
         ("validateImpossToCalc_mThreshold", int, 3),
         ("in:SyntenyBlocks", str, 'None'),
         ("out:SyntenyBlocks", str, "./res/syntenyBlocksDrawer.txt"),
         ("mode:chromosomesRewrittenInTbs", bool, False),
         ('convertGenicToTbCoordinates', bool, False),
         ("distanceMetric", str, 'CD'),
         ('nbHpsRecommendedGap', int, 2),
         ('targetProbaRecommendedGap', float, 0.01),
         ("scaleFactorRectangles", float, 2.0),
         ("out:ImageName", str, "./res/homologyMatrix.svg"),
         ('verbose', bool, True)],
        __doc__)

    # Load genomes
    genome1 = myLightGenomes.LightGenome(arguments["genome1"], withDict=True)
    genome2 = myLightGenomes.LightGenome(arguments["genome2"], withDict=True)
    # Change genome format
    genome1Name = genome1.name
    genome2Name = genome2.name
    sbsInPairComp = None
    if arguments['in:SyntenyBlocks'] is not 'None':
        # load synteny blocks
        sbsInPairComp = myDiags.parseSbsFile(arguments['in:SyntenyBlocks'], genome1=genome1, genome2=genome2)
    families = myLightGenomes.Families(arguments["families"])
    modesFilter = list(myDiags.FilterType._keys)
    filterType = myDiags.FilterType[modesFilter.index(arguments["filterType"])]

    tandemGapMax = arguments['tandemGapMax']
    gapMax = arguments['gapMax']
    identifyBreakpointsWithinGaps = arguments['identifyBreakpointsWithinGaps']
    nonOverlappingSbs = arguments['nonOverlappingSbs']
    overlapMax = arguments['overlapMax']
    consistentSwDType = arguments['consistentSwDType']
    considerMonogenicSb=arguments['considerMonogenicSb']
    minChromLength = arguments['minChromLength']
    pThreshold = arguments['pThreshold']
    validateImpossToCalc_mThreshold = arguments['validateImpossToCalc_mThreshold']
    chromosomesRewrittenInTbs = arguments['mode:chromosomesRewrittenInTbs']
    convertGenicToTbCoordinates = arguments['convertGenicToTbCoordinates']
    distanceMetric = arguments['distanceMetric']
    nbHpsRecommendedGap = arguments['nbHpsRecommendedGap']
    targetProbaRecommendedGap = arguments['targetProbaRecommendedGap']
    scaleFactorRectangles = arguments['scaleFactorRectangles']

    #thresholdChr = 50

    #by convention:
    if gapMax == 'None':
        gapMax = None
    else:
        try:
            gapMax = int(gapMax)
        except:
            raise ValueError('gapMax must be an int or None')

    assert distanceMetric == 'DPD' or distanceMetric == 'MD'\
        or distanceMetric == 'CD' or distanceMetric == 'ED'
    assert (convertGenicToTbCoordinates and chromosomesRewrittenInTbs)\
        or not convertGenicToTbCoordinates

    if not chromosomesRewrittenInTbs:
        (chr1, range1) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr1:deb1-fin1"], genome1)
        (chr2, range2) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr2:deb2-fin2"], genome2)
        chrom1 = myLightGenomes.LightGenome()
        chrom2 = myLightGenomes.LightGenome()
        chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
        chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]

        if sbsInPairComp is not None:
            new_sbsInPairComp = myTools.Dict2d(list)
            for sb in sbsInPairComp[chr1][chr2]:
                newl1 = []
                newl2 = []
                newla = []
                for (idxHp, aG) in enumerate(sb.la):
                    tb1 = []
                    tb2 = []
                    for i1g in sb.l1[idxHp]:
                        if (range1[0] <= i1g and i1g <= range1[1]):
                            tb1.append(i1g)
                    for i2g in sb.l2[idxHp]:
                        if (range2[0] <= i2g and i2g <= range2[1]):
                            tb2.append(i2g)
                    if len(tb1) > 0 and len(tb2) > 0:
                        newl1.append(tb1)
                        newl2.append(tb2)
                        newla.append(aG)
                new_sbsInPairComp[chr1][chr2].append(myDiags.SyntenyBlock(myDiags.Diagonal(sb.dt, newl1, newl2, newla), sb.pVal))
            sbsInPairComp = new_sbsInPairComp
        else:
            # extract diagonals in the ROI without considering other pairwise
            # comparisons
            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(chrom1,
                                                                chrom2,
                                                                families,
                                                                tandemGapMax=tandemGapMax,
                                                                gapMax=gapMax,
                                                                identifyBreakpointsWithinGaps=identifyBreakpointsWithinGaps,
                                                                nonOverlappingSbs=nonOverlappingSbs,
                                                                overlapMax=overlapMax,
                                                                consistentSwDType=consistentSwDType,
                                                                filterType=filterType,
                                                                considerMonogenicSb=considerMonogenicSb,
                                                                minChromLength=minChromLength,
                                                                distanceMetric=distanceMetric,
                                                                pThreshold=pThreshold,
                                                                validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold)

        genesDiagIndices = []
        for sb in sbsInPairComp[chr1][chr2]:
            genesDiagIndices.append([])
            for idxHp, aG in enumerate(sb.la):
                for (gi1, gi2) in itertools.product(sb.l1[idxHp], sb.l2[idxHp]):
                    genesDiagIndices[-1].append((gi1, gi2))

        ###
        # Build Genes Strands
        ###
        genesStrandsC1 = [s for (_, s) in chrom1[chr1]]
        genesStrandsC2 = [s for (_, s) in chrom2[chr2]]

        ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
         genesHomologiesHpSign,
         (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
         genesHomologyGroupsInWindow) =\
            drawHomologyMatrixWithSBs.genesComputeHomologyInformations(chr1, chr2, chrom1, chrom2,
                                                                       families,
                                                                       filterType,
                                                                       minChromLength,
                                                                       tandemGapMax)


        strArray = drawHomologyMatrixWithSBs.drawHomologyMatrix((range1, range2),
                                                                (genesStrandsC1, genesStrandsC2),
                                                                (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
                                                                (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
                                                                genesHomologiesHpSign,
                                                                genesHomologyGroupsInWindow,
                                                                genesDiagIndices,
                                                                outputFileName=arguments["out:ImageName"],
                                                                maxWidth=100,
                                                                maxHeight=100,
                                                                scaleFactorRectangles=scaleFactorRectangles)

    else:
        #chromosomes are shown as lists of tbs
        g1fId = myMapping.labelWithFamID(genome1, families)
        g2fId = myMapping.labelWithFamID(genome2, families)
        ((g1f, gf2gfId1, (nCL1, nGL1)),
         (g2f, gf2gfId2, (nCL2, nGL2))) =\
            myDiags.filter2D(g1fId, g2fId,
                             filterType,
                             minChromLength,
                             keepOriginal=True)
        (g1tb, gtb2gf1, nGTD1) =\
            myMapping.remapRewriteInTb(g1f,
                                       tandemGapMax=tandemGapMax,
                                       mOld=None)
        (g2tb, gtb2gf2, nGTD2) =\
            myMapping.remapRewriteInTb(g2f,
                                       tandemGapMax=tandemGapMax,
                                       mOld=None)
        gtb2gfId1 = {}
        gfId2gtb1 = {}
        for c in gtb2gf1:
            # see Mapping class addition
            gtb2gfId1[c] = gtb2gf1[c] + gf2gfId1[c]
            gfId2gtb1[c] = gtb2gfId1[c].old
        gtb2gfId2 = {}
        gfId2gtb2 = {}
        for c in gtb2gf2:
            # see Mapping class addition
            gtb2gfId2[c] = gtb2gf2[c] + gf2gfId2[c]
            gfId2gtb2[c] = gtb2gfId2[c].old

        if not convertGenicToTbCoordinates:
            (chr1, range1) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr1:deb1-fin1"], g1tb)
            (chr2, range2) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr2:deb2-fin2"], g2tb)
        else:
            (chr1, range1) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr1:deb1-fin1"], g1tb, g2gtb=gfId2gtb1)
            (chr2, range2) = drawHomologyMatrixWithSBs.parseChrRange(arguments["chr2:deb2-fin2"], g2tb, g2gtb=gfId2gtb2)

        chrom1_tb = {}
        chrom2_tb = {}
        chrom1_tb[chr1] = g1tb[chr1][range1[0]:range1[1]]
        chrom2_tb[chr2] = g2tb[chr2][range2[0]:range2[1]]

        #Focus on the chromosome of the window, just give simple name to the chromosome of interest
        chrom1 = genome1[chr1]
        Ctb2Cf1 = gtb2gf1[chr1]
        CfId2Ctb1 = gfId2gtb1[chr1]

        chrom2 = genome2[chr2]
        Ctb2Cf2 = gtb2gf2[chr2]
        CfId2Ctb2 = gfId2gtb2[chr2]

        ###
        # Build TbNumberOfGenesInEachTbC1 : [ 4,5,1,1,6,2, ...] number og genes in each TB of C1
        ###
        TbNumberOfGenesInEachTbC1 = [len(Ctb2Cf1[i1_tb]) for i1_tb in range(range1[0], range1[1])]
        TbNumberOfGenesInEachTbC2 = [len(Ctb2Cf2[i2_tb]) for i2_tb in range(range2[0], range2[1])]

        ###
        # Build TBStrands
        ###
        TbStrandsC1 = [s for (_, s) in chrom1_tb[chr1]]
        TbStrandsC2 = [s for (_, s) in chrom2_tb[chr2]]

        ###
        # Build rangeXTB
        ###
        (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2),
         TbHomologyGroupsInWindow) =\
            drawHomologyMatrixWithSBs.TbComputeHomologyInformations(chrom1_tb[chr1], chrom2_tb[chr2])
        ###
        # Convert into the correct format for the function TbHomologyGroupsInWindow
        ###
        tmpTbHomologyGroupsInWindow = []
        for (tbs1, tbs2) in TbHomologyGroupsInWindow:
            tmpTbs1 = []
            tmpTbs2 = []
            for tb1 in tbs1:
                tmpTbs1.append([tb1])
            for tb2 in tbs2:
                tmpTbs2.append([tb2])
            tmpTbHomologyGroupsInWindow.append((tmpTbs1, tmpTbs2))
        TbHomologyGroupsInWindow = tmpTbHomologyGroupsInWindow
        TbNoHomologiesInWindowC1 = [[tb1] for tb1 in TbNoHomologiesInWindowC1]
        TbNoHomologiesInWindowC2 = [[tb2] for tb2 in TbNoHomologiesInWindowC2]

        if sbsInPairComp is not None:
            new_sbsInPairComp = myTools.Dict2d(list)
            for sb in sbsInPairComp[chr1][chr2]:
                # change the sb.lX structure from list of lists to list of ints
                new_l1 = [CfId2Ctb1[tb[0]] for tb in sb.l1]
                new_l2 = [CfId2Ctb2[tb[0]] for tb in sb.l2]
                sb = myDiags.SyntenyBlock(myDiags.Diagonal(sb.dt, new_l1, new_l2, sb.la), sb.pVal)
                if (range1[0] <= sb.minOnG(1) and sb.maxOnG(1) <= range1[1])\
                        and (range2[0] <= sb.minOnG(2) and sb.maxOnG(2) <= range2[1]):
                    # sb is perfectly included in the ROI
                    new_sbsInPairComp[chr1][chr2].append(sb)
                elif (sb.maxOnG(1) < range1[0] or range1[1] < sb.minOnG(1))\
                        or (sb.maxOnG(2) < range2[0] or range2[1] < sb.minOnG(2)):
                    # sb is not in the ROI
                    continue
                else:
                    # sb is partially included in the ROI
                    sb.truncate(range1, range2)
                    new_sbsInPairComp[chr1][chr2].append(sb)
            sbsInPairComp = new_sbsInPairComp
        else:
            # extract diagonals in the ROI without considering other pairwise
            # comparisons
            sbsInPairComp = myDiags.extractSbsInPairCompGenomesInTbs(chrom1_tb, chrom2_tb,
                                                                     gapMax=gapMax,
                                                                     distanceMetric=distanceMetric,
                                                                     pThreshold=pThreshold,
                                                                     identifyBreakpointsWithinGaps=identifyBreakpointsWithinGaps,
                                                                     nonOverlappingSbs=nonOverlappingSbs,
                                                                     overlapMax=overlapMax,
                                                                     consistentSwDType=consistentSwDType,
                                                                     considerMonogenicSb=considerMonogenicSb,
                                                                     nbHpsRecommendedGap=nbHpsRecommendedGap,
                                                                     targetProbaRecommendedGap=targetProbaRecommendedGap,
                                                                     validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold)

        TbDiagIndices = []
        for sb in sbsInPairComp[chr1][chr2]:
            TbDiagIndices.append([])
            for idxHp, aG in enumerate(sb.la):
                assert isinstance(sb.l1[idxHp], int) and isinstance(sb.l2[idxHp], int)
                TbDiagIndices[-1].append((sb.l1[idxHp], sb.l2[idxHp]))

        strArray =\
            drawHomologyMatrixWithSBs.drawHomologyMatrix((range1, range2),
                                                         (TbStrandsC1, TbStrandsC2),
                                                         ([], []),
                                                         (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2),
                                                         TbHpSign,
                                                         TbHomologyGroupsInWindow,
                                                         TbDiagIndices,
                                                         outputFileName=arguments["out:ImageName"],
                                                         maxWidth=100,
                                                         maxHeight=100,
                                                         symbolsInGenes=(TbNumberOfGenesInEachTbC1, TbNumberOfGenesInEachTbC2),
                                                         scaleFactorRectangles=scaleFactorRectangles
                                                         )

    #copy the css style sheet
    dirNameImage = os.path.dirname(arguments["out:ImageName"])
    dirNameImage = dirNameImage if dirNameImage != "" else "."
    print >> sys.stderr, "cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage)
    os.system("cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage))

    # write a simple file with all diagonals into output file
    f = open(arguments['out:SyntenyBlocks'], 'w')
    print >> f, "Mode : %s" % 'Genic scale' if arguments['mode:chromosomesRewrittenInTbs'] is False else 'Tandem Blocks scale'
    print >> f, "chromosome %s de %s\t%s\t%s\tchromosome %s de %s\t%s\t%s\t%s" % (chr1, genome1Name, 'beginC1', 'endC1', chr2, genome2Name, 'beginC2', 'endC2', 'length in families')
    print >> f, "c1\tbeg1\tend1\tc2\tbeg2\tend2\thps\tpVal"

    for sb in sbsInPairComp[chr1][chr2]:
        if isinstance(sb.l1[0], list) and isinstance(sb.l2[0], list):
            minl1 = min(sb.l1[0])
            maxl1 = max(sb.l1[-1])
            minl2 = min(sb.l2[0])
            maxl2 = max(sb.l2[-1])
        elif isinstance(sb.l2[0], int) and isinstance(sb.l2[0], int):
            minl1 = min(sb.l1)
            maxl1 = max(sb.l1)
            minl2 = min(sb.l2)
            maxl2 = max(sb.l2)
        # indices of genes start at 1 to be coherent with the output image
        print >> f, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chr1, minl1 + 1, maxl1 + 1, chr2, minl2 + 1, maxl2 + 1, len(sb.la), sb.pVal)
    f.close()

    # Add lengends and title to the ouput matrix
    height = 100
    width = 100
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
    if not arguments['mode:chromosomesRewrittenInTbs']:
        title =\
            "%s, tandemGapMax=%s tbs, gapMax=%s%s, overlapMax=%s tbs, %s sbs" %\
            ('MH',
             arguments['tandemGapMax'],
             arguments['gapMax'],
             arguments['distanceMetric'],
             arguments['overlapMax'],
             len(list(sbsInPairComp.iteritems2d())))
    else:
        title =\
            "%s, tandemGapMax=%s tbs, gapMax=%s%s, overlapMax=%s tbs, %s sbs" %\
            ('MHP',
             arguments['tandemGapMax'],
             arguments['gapMax'],
             arguments['distanceMetric'],
             arguments['overlapMax'],
             len(list(sbsInPairComp.iteritems2d())))\

    var += ['<svg x="5" y="0" viewBox="5 0 95 5" width="95" height="5" xmlns="http://www.w3.org/2000/svg" xmlns:xhtml="http://www.w3.org/1999/xhtml" xmlns:xlink="http://www.w3.org/1999/xlink">\n',
                    '<foreignObject x="0" y="0" width="95" height="5">\n',
                    '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
                            '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
                            '<xhtml:div style="color:black; word-wrap:break-word; font-size:1.5px; font-family:Arial" >' + title + '\n',
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

#Find diags with the more paralogs
@myTools.deprecated
def searchInterestingDiags(listOfDiags, range1, range2):
    ecart=-sys.maxint-1
    diag_=None
    for diag in listOfDiags:
        l1=diag[0][1]
        l2=diag[1][1]
        la=diag[2]
        if len(l1) != len(l2):
            print >> sys.stderr, "assymetric (on genome1 and genome2) SB of ", len([ anc for (anc,_,_,_) in la]), " HPs"
            print >> sys.stderr, "diag on G1 from %s to %s" % (min(l1[0][0]+range1[0], l1[-1][0]+range1[0]) , max(l1[0][0]+range1[0], l1[-1][0]+range1[0]))
            print >> sys.stderr, "diag on G2 from %s to %s" % (min(l2[0][0]+range2[0], l2[-1][0]+range2[0]) , max(l2[0][0]+range2[0], l2[-1][0]+range2[0]))
            if ecart < abs(len(l1) -  len(la)):
                ecart = abs(len(l1) -  len(la))
                diag_=diag
            if ecart < abs(len(l2) -  len(la)):
                ecart = abs(len(l2) -  len(la))
                diag_=diag
    diag = diag_ if diag_ != None else listOfDiags[0]
    l1 = diag[0][1]
    l2 = diag[1][1]
    la = diag[2]
    print >> sys.stderr, "The most assymetric diag"
    print >> sys.stderr, " SB of ", len([ anc for (anc,_,_,_) in la]), " HPs"
    print >> sys.stderr, "diag on G1 = ", [ gene[0]+range1[0] for gene in l1]
    print >> sys.stderr, "diag on G2 = ", [ gene[0]+range2[0] for gene in l2]

    minIndiceG1 = sys.maxint
    maxIndiceG1 = -sys.maxint-1
    minIndiceG2 = sys.maxint
    maxIndiceG2 = -sys.maxint-1
    for diag in listOfDiags :
        ((_,l1),(_,l2),_) = diag
        for (i1,_) in l1:
            minIndiceG1 = i1 if i1 < minIndiceG1 else minIndiceG1
            maxIndiceG1 = i1 if i1 > maxIndiceG1 else maxIndiceG1
        for (i2,_) in l2:
            minIndiceG2 = i2 if i2 < minIndiceG2 else minIndiceG2
            maxIndiceG2 = i2 if i2 > maxIndiceG2 else maxIndiceG2
    print >> sys.stderr, 'Indices entre lesquels se trouvents les diagonales sur le G1 = [%s,%s]' % (minIndiceG1+range1[0], maxIndiceG1+range1[0])
    print >> sys.stderr, 'Indices entre lesquels se trouvents les diagonales sur le G2 = [%s,%s]' % (minIndiceG2+range2[0], maxIndiceG2+range2[0])

    return

# ask the user for the desired chromosome ranges
@myTools.deprecated
def chooseChrsAndRanges(genome1, genome2, families, distanceMetric = 'DPD'):
    while True:
        try:
            (chr1, range1) = drawHomologyMatrixWithSBs.parseChrRange(raw_input("chr1:deb1-fin1 = "), genome1)
            break
        except ValueError:
            print >> sys.stderr, "You need to write somtehing as chr1:deb1-fin1 with chr1 a chr of G1 and deb1 and fin1 indices of the first and last gene (indices start at 1)"

    while True:
        try:
            (chr2, range2) = drawHomologyMatrixWithSBs.parseChrRange(raw_input("chr2:deb2-fin2 = "), genome2)
            break
        except ValueError:
            print >> sys.stderr, "You need to write somtehing as chr2:deb2-fin2 with chr2 a chr of G2 and deb2 and fin2 indices of the first and last gene (indices start at 1)"

    chrom1 ={}
    chrom2 ={}
    chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
    chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]
    listOfDiags = myDiags.extractSbsInPairCompGenomes(chrom1, chrom2, families, gapMax=arguments["gapMax"], consistentSwDType=arguments["consistentSwDType"], filterType=filterType, minChromLength=arguments["minChromLength"], distanceMetric=arguments["distanceMetric"], verbose=arguments['verbose'])
    listOfDiags = list(listOfDiags)
    print >> sys.stderr, "pairwise comparison of the two chromosomes yields" , len(listOfDiags), "diagonals."
    if len(listOfDiags) == 0:
        print >> sys.stderr, "There is no diag in the considered region of interest"
    else:
        #recherche des diagonales interessantes a afficher
        print >> sys.stderr, range1
        print >> sys.stderr, range2
        searchInterestingDiags(listOfDiags, range1, range2)

    print >> sys.stderr, "Do you want to chose a new region of interest ?"
    if raw_input('y/n ? ') == 'y':
        return chooseChrsAndRanges(genome1, genome2, families)
    else:
        return (chrom1, chrom2, range1, range2)
