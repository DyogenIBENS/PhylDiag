#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

"""
Build the homology matrix with synteny blocks using the mySvgDrawer Library
"""

import sys
import collections
import itertools
import random
import utils.myTools as myTools
import utils.mySvgDrawer as svgDrw
import utils.myDiags as myDiags
import utils.myMapping as myMapping


# parse the user input (text) for the chromosome range and asses if this query is consistent with the genome data
# g2gtb is a dictionnary to convert indices from gene coordinates into tb coordinates
def parseChrRange(text, genome, g2gtb=None):
    if len(text.split(":")) == 2:
        chr = int(text.split(":")[0]) if text.split(":")[0].isdigit() else text.split(":")[0]
        if chr not in genome.keys():
            print >> sys.stderr, "chr %s not in genome" % chr
            raise ValueError
    else:
        print >> sys.stderr, "range not formated as expected : \"chr:deb-fin\""
        raise ValueError

    range = text.split(":")[-1].split("-")
    if len(range) == 2 and range[0].isdigit and (range[1].isdigit() or range[1] == '~'):
        if g2gtb is not None:
            range[0] = g2gtb[chr][int(range[0])-1]
            range[1] = g2gtb[chr][int(range[1])-1] if range[1] != '~' else len(genome[chr])
        else:
            range = (int(range[0])-1, int(range[1]) if range[1] != '~' else len(genome[chr]))
        if range[0] < 0 or range[1] > len(genome[chr]) or range[1] <= 0 or range[1] <= range[0]:
            print >> sys.stderr,\
                "range %s is incoherent for chr %s. FYI chr %s contains %s elements. Be sure that beginning < end of range" % ([range[0]+1, range[1]], chr, chr, len(genome[chr]))
            raise ValueError
    else:
        raise ValueError
    return (chr, range)


def TbComputeHomologyInformations(chrom1_tb, chrom2_tb):
    ###
    # Build MHP = { i1_tb : {i2_tb : hpSign, ...}...} with hpSign
    # the hpSign (s1*s2) of the homology at the (i1,i2) coordinate
    ###
    (MHP, locG2) = myDiags.homologyMatrix(chrom1_tb, chrom2_tb)

    tbsOnChr2ThatHaveHomologies = set([])
    for tb1 in MHP:
        tbsOnChr2ThatHaveHomologies |= set(MHP[tb1].keys())

    ###
    # Build TBNoHomologiesInWindowX = [..., i_tb, ... ] with i_tb the index of
    # the tb of chromosomeX_tb with no homology in the window
    ###
    TbNoHomologiesInWindowC1 = []
    for (i1_tb, _) in enumerate(chrom1_tb):
        if i1_tb not in MHP:
            TbNoHomologiesInWindowC1.append(i1_tb)
    TbNoHomologiesInWindowC2 = []
    for (i2_tb, _) in enumerate(chrom2_tb):
        if i2_tb not in tbsOnChr2ThatHaveHomologies:
            TbNoHomologiesInWindowC2.append(i2_tb)
    ###
    # Build homologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tbCX_X : [..., i, ...] list of indices of genes in the tb
    ###
    locC1_tbIndx = collections.defaultdict(list)
    locC2_tbIndx = collections.defaultdict(list)
    TbHomologyGroupsInWindow = []
    listOfHPsCoordinates = []
    for i1_tb in MHP:
        for i2_tb in MHP[i1_tb]:
            listOfHPsCoordinates.append((i1_tb, i2_tb))
    for (i1_tb, i2_tb) in listOfHPsCoordinates:
        locC2_tbIndx[i1_tb].append(i2_tb)
        locC1_tbIndx[i2_tb].append(i1_tb)
    HomologyGroup1 = set([])
    HomologyGroup2 = set([])
    for (i1_tb, i2_tb) in listOfHPsCoordinates:
        if i1_tb not in HomologyGroup1 or i2_tb not in HomologyGroup2:
            TbHomologyGroupsInWindow.append(([], []))
            HomologyGroup1.add(i1_tb)
            for i2_tbb in locC2_tbIndx[i1_tb]:
                TbHomologyGroupsInWindow[-1][1].append(i2_tbb)
                HomologyGroup2.add(i2_tbb)
                for i1_tbb in locC1_tbIndx[i2_tbb]:
                    TbHomologyGroupsInWindow[-1][0].append(i1_tbb)
                    HomologyGroup1.add(i1_tbb)
                del locC1_tbIndx[i2_tbb]
            del locC2_tbIndx[i1_tb]
    del HomologyGroup1
    del HomologyGroup2

    return (MHP, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow)


def genesComputeHomologyInformations(chr1, chr2, chrom1, chrom2, ancGenes,
                                     filterType,
                                     minChromLength, tandemGapMax):
    c1_aID = myMapping.labelWithAncGeneID(chrom1, ancGenes)
    c2_aID = myMapping.labelWithAncGeneID(chrom2, ancGenes)
    # Must be applied on the two genomes
    ((c1_aID_filt, Cf2CaID1, (nCL1, nGL1)),
     (c2_aID_filt, Cf2CaID2, (nCL2, nGL2))) = myDiags.filter2D(c1_aID, c2_aID,
                                                               filterType,
                                                               minChromLength,
                                                               keepOriginal=True)
    (chrom1_tb, Ctb2Cf1, nGTD1) = myMapping.remapRewriteInTb(c1_aID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=None)
    (chrom2_tb, Ctb2Cf2, nGTD2) = myMapping.remapRewriteInTb(c2_aID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=None)
    Ctb2CaID1 = {}
    for c in Ctb2Cf1:
        # see Mapping class addition
        Ctb2CaID1[c] = Ctb2Cf1[c] + Cf2CaID1[c]
    Ctb2CaID2 = {}
    for c in Ctb2Cf2:
        # see Mapping class addition
        Ctb2CaID2[c] = Ctb2Cf2[c] + Cf2CaID2[c]

    #Focus on the chromosome of the window
    chrom1_ = chrom1[chr1]
    chrom2_ = chrom2[chr2]
    c1_aID = c1_aID[chr1]
    c2_aID = c2_aID[chr2]
    c1_aID_filt = c1_aID_filt[chr1]
    c2_aID_filt = c2_aID_filt[chr2]
    Cf2CaID1 = Cf2CaID1[chr1]
    Cf2CaID2 = Cf2CaID2[chr2]
    chrom1_tb = chrom1_tb[chr1]
    chrom2_tb = chrom2_tb[chr2]
    Ctb2Cf1 = Ctb2Cf1[chr1]
    Ctb2Cf2 = Ctb2Cf2[chr2]
    Ctb2CaID1 = Ctb2CaID1[chr1]
    Ctb2CaID2 = Ctb2CaID2[chr2]

    ###
    # Build genesRemovedDuringFilteringCX = [..., i, ...] the list of genes that
    # have been removed during the filtering process
    ###
    genesRemovedDuringFilteringC1 = [i1 for (i1, (anc, _)) in enumerate(c1_aID) if anc is None]
    genesRemovedDuringFilteringC2 = [i2 for (i2, (anc, _)) in enumerate(c2_aID) if anc is None]

    (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow) =\
        TbComputeHomologyInformations(chrom1_tb, chrom2_tb)

    genesHomologiesHpSign = collections.defaultdict(lambda: collections.defaultdict(int))
    for i1_tb in TbHpSign:
        for i2_tb in TbHpSign[i1_tb]:
            for (i1, i2) in itertools.product([ii1 for ii1 in Ctb2CaID1[i1_tb]],
                                              [ii2 for ii2 in Ctb2CaID2[i2_tb]]):
                s1 = chrom1_[i1][1]
                s2 = chrom2_[i2][1]
                genesHomologiesHpSign[i1][i2] = s1*s2

    ###
    # Build genesNoHomologiesInWindow1 = [..., [i5,i6,i7], ...] list of tbs with [i5,i6,i7] a tb of three genes whose indices are i5,i6 and i7
    ###
    genesNoHomologiesInWindowC1 = []
    genesNoHomologiesInWindowC2 = []
    for i1_tb in TbNoHomologiesInWindowC1:
        genesNoHomologiesInWindowC1.append([i1 for i1 in Ctb2CaID1[i1_tb]])
    for i2_tb in TbNoHomologiesInWindowC2:
        genesNoHomologiesInWindowC2.append([i2 for i2 in Ctb2CaID2[i2_tb]])

    ###
    # Build genesHomologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tb_CX_X : [..., i, ...] list of indices of genes in the tb
    ###
    genesHomologyGroupsInWindow = []
    for (tbs1, tbs2) in TbHomologyGroupsInWindow:
        genesHomologyGroupsInWindow.append(([], []))
        for i1_tb_group in tbs1:
            genesHomologyGroupsInWindow[-1][0].append([i1 for i1 in Ctb2CaID1[i1_tb_group]])
        for i2_tb_group in tbs2:
            genesHomologyGroupsInWindow[-1][1].append([i2 for i2 in Ctb2CaID2[i2_tb_group]])

    return ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
            genesHomologiesHpSign,
            (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
            genesHomologyGroupsInWindow)


# Generator of levels for colors or gray indices within a palette:
class levelIdxGenerator():
    def __init__(self, farIdxs=False, grays=False):
        if grays is not None:
            # see the css joint file "HomologyGroup" ranges from 0 to 44 (included)
            self.firstLevelIdx = 2
            self.lastLevelIdx = 44
        else:
            # see the css joint file "NoHomologyInWindow" ranges from 0 to 14 (included)
            self.firstLevelIdx = 2
            self.lastLevelIdx = 14
        if farIdxs is False:
            self.availableLevels = range(self.firstLevelIdx, self.lastLevelIdx+1)
        else:
            self.availableLevels = range(self.firstLevelIdx, self.lastLevelIdx+1, 3) +\
                                   range(self.firstLevelIdx+1, self.lastLevelIdx+1, 3) +\
                                   range(self.firstLevelIdx+2, self.lastLevelIdx+1, 3)
        self.currIdx = 0
        assert len(self.availableLevels) == (self.lastLevelIdx+1) - self.firstLevelIdx

    def getLevel(self, differentFrom=set([])):
        if len(differentFrom.intersection(self.availableLevels)) == len(self.availableLevels):
            print >> sys.stderr, "Warning: too many colors too avoid, thus a non-optimal choice of the color is made"
        else:
            while self.currIdx in differentFrom:
                if self.currIdx < self.lastLevelIdx - self.firstLevelIdx:
                    self.currIdx += 1
                else:
                    self.currIdx = 0
        level = self.availableLevels[self.currIdx]
        self.currIdx += 1
        if self.currIdx > self.lastLevelIdx - self.firstLevelIdx:
            self.currIdx = 0
        return level


def neighboursLevels(chromosome, i):
    # level of the 'L'eft neighbour
    levelL = chromosome[i-1].SVGclass if i-1 in chromosome else None
    # level of the 'R'ight neighbour
    levelR = chromosome[i+1].SVGclass if i+1 in chromosome else None

    def convertSVGclassIntoInt(strSVGclass):
        if strSVGclass is not None:
            try:
                int(strSVGclass[-2:])
            except:
                int(strSVGclass[-1])

    levelL = convertSVGclassIntoInt(levelL)
    levelR = convertSVGclassIntoInt(levelR)
    neighboursLevels = set([levelL, levelR])
    neighboursLevels = set([l for l in neighboursLevels if l is not None])
    return neighboursLevels


# draw either the mh or the mhp, if draw mode is 'writeinTB'
# inputs :
#       genesStrandsCX = [+1, -1, ...] of length = to nX
#       genesRemovedDuringFilteringC1 = [..., i, ...] with i the index of the gene removed during the filtering process (CX : Chromosome X)
#      tbWithNoHomologyInWindowC1 = [..., [i6,i7,i8], ...] list of tbs with no homologies in the window, inside the index of genes. If draw mode is 'writeinTB' : tbWithNoHomologyInWindowC1 = [..., i6, ...] : just the index of the TB
#       hpSigns = { i1 : {i2 : s, ...}...} with hpSign the hp sign (s1*s2) of the homology at the (i1,i2) coordinate (i1-th gene on C1 and i2-th gene on C2)
#      homologyGroupsInWindow = [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...] a list of 2-uples of homologous tbs indices to be coloured in the same color
#              with tb_CX_X : [..., i, ...] list of indices of genes in the tb
#       diagIndices = [..., [...,(i16,j16),...], ...] list of diagonals with diagonals = list of all the points of the diagonal
# output :
#       string with the svg drawing of the mhp (or the mh)
def drawHomologyMatrix(((begC1, endC1), (begC2, endC2)), (genesStrandsC1, genesStrandsC2),
                       (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
                       (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2),
                       hpSigns, homologyGroupsInWindow, diagsIndices,
                       outputFileName=None, maxWidth=100, maxHeight=100, symbolsInGenes=None):
    # For the print on screen
    begC1 = begC1 + 1
    begC2 = begC2 + 1
    endC1 = endC1 + 1
    endC2 = endC2 + 1
    # nx = number of genes on the first genome
    # ny = number of genes on the second genome
    nx = len(genesStrandsC1)
    ny = len(genesStrandsC2)

    # the size of the components of the matrix is chosen using the smallest and more restricting dimension (contains more genes comparing to its size)
    sizeCase = float(min(float(maxWidth) / (nx + 3), float(maxHeight) / (ny + 3)))  # +3 for margins
    sizeText = float(sizeCase*0.9)
    width = float(nx+3) * sizeCase
    height = float(ny+3) * sizeCase
    print >> sys.stderr, "width=", width
    print >> sys.stderr, "height=", height

    scene = svgDrw.Scene(name='homology_matrix', width=width, height=height)

    nbLinesX = nx+1 #Nb of vertical lines (x varies) in the matrix
    nbLinesY = ny+1 #Nb of horizontal lines (y varies) in the matrix

    farColorsGenerator = levelIdxGenerator(farIdxs=True)
    closeColorsGenerator = levelIdxGenerator()
    farGraysGenerator = levelIdxGenerator(farIdxs=True, grays=True)

    # draw lines of chromosomes
    # offset_genes : corresponds to the chromosome lines positions passing through the middle of the genes
    offset_genes_x = sizeCase + sizeCase/2
    offset_genes_y = sizeCase + sizeCase/2
    scene.add(svgDrw.Line((offset_genes_x, height - offset_genes_y), (width-sizeCase, height - offset_genes_y), width=0.1*sizeCase))
    scene.add(svgDrw.Line((offset_genes_x, sizeCase), (offset_genes_x, height - (offset_genes_y)), width=0.1*sizeCase))
    offset_matrix_x = 2 * sizeCase
    offset_matrix_y = 2 * sizeCase
    # draw Diagonals first because they are on the background
    print >> sys.stderr, "Nb of diagonals showed = ", len(diagsIndices)
    # sort diagonals to give colors according to localisations of diagonals
    diagsIndices.sort(key=lambda x: x[0])
    for diag in diagsIndices:
        # choose a color different from the neighbours
        color = farColorsGenerator.getLevel()
        for (i, j) in diag:
            cx_s = i*sizeCase
            cy_s = j*sizeCase
            scene.add(svgDrw.Rectangle((cx_s+offset_matrix_x, height-(cy_s+sizeCase+offset_matrix_y)), sizeCase, sizeCase, fill_opacity=0.90, svgClass="HomologGroup%s" % color))
            #scene.add(scgDrw.Rectangle((cx_s,height-(cy_s+sizeCase)), sizeCase, sizeCase, fill_opacity=0.90, svgClass = "chromgrp%s" % color ))
            #scene.add(svgDrw.Rectangle((cx_s,height-(cy_s+sizeCase)), sizeCase, sizeCase, fill_opacity=fill_opacity=0.90, color=(color,color,color))
        # draw rectangles around diagonals
        min_i = min(diag, key=lambda x: x[0])[0]
        max_i = max(diag, key=lambda x: x[0])[0]
        min_j = min(diag, key=lambda x: x[1])[1]
        max_j = max(diag, key=lambda x: x[1])[1]
        cx_s = min_i*sizeCase
        cy_s = max_j*sizeCase
        scene.add(svgDrw.Rectangle((cx_s+offset_matrix_x, height-(cy_s+offset_matrix_y+sizeCase)), (max_j-min_j)*sizeCase + sizeCase, (max_i-min_i)*sizeCase + sizeCase,  stroke='black', fill='none', strokeWidth=0.2*sizeCase))

    # tick lines
    widthTicks = min(float(width)/1000, float(height)/1000)
    sizeTextTicks = widthTicks*10
    for (i, ni) in enumerate(range(begC1, endC1)):
        cx = i*sizeCase
        if ni % 10 == 0:
            scene.add(svgDrw.Line((offset_matrix_x+sizeCase/2+cx, height - offset_genes_y/2), (offset_matrix_x+sizeCase/2+cx, height - offset_genes_y), width=widthTicks))
        if ni % 50 == 0:
            cxText = offset_matrix_x+sizeCase/2+cx
            cyText = height - max(offset_genes_y/2, sizeTextTicks/2)
            cxx = offset_matrix_x+cx+sizeCase/2
            if nx > 750 or ny > 750:
                if ni % 100 == 0:
                    scene.add(svgDrw.Text((cxText, cyText), str(ni), text_anchor="middle", size=sizeTextTicks))
                    scene.add(svgDrw.Line((cxx, height - offset_matrix_y), (cxx, sizeCase), width=sizeCase*0.1))
            else:
                scene.add(svgDrw.Text((cxText, cyText), str(ni), text_anchor="middle", size=sizeTextTicks))
                scene.add(svgDrw.Line((cxx, height - offset_matrix_y), (cxx, sizeCase), width=sizeCase*0.1))

    for (j, nj) in enumerate(range(begC2, endC2)):
        cy = j*sizeCase
        if nj % 10 == 0:
            scene.add(svgDrw.Line((offset_genes_x/2, height - (offset_matrix_y+sizeCase/2+cy)), (offset_genes_x, (height - (offset_matrix_y+sizeCase/2+cy))), width=widthTicks))
        if nj % 50 == 0:
            cyText = height - (offset_matrix_y+sizeCase/2+cy)
            cxText = max(offset_genes_x/2, sizeTextTicks/2)
            cxx = offset_genes_x/2
            cyy = height - (offset_matrix_y+sizeCase/2+cy)
            if nx > 750 or ny > 750:
                if nj % 100 == 0:
                    #scene.add(svgDrw.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCase,cxx+sizeCase,cyy)))
                    scene.add(svgDrw.Text((cxText, cyText), str(nj), text_anchor="middle", size=sizeTextTicks, transform="rotate(90,%s,%s)" % (cxText, cyText)))
                    scene.add(svgDrw.Line((offset_matrix_x, cyy), (width-sizeCase, cyy), width=sizeCase*0.1))
            else:
                #scene.add(svgDrw.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCase,cxx+sizeCase,cyy)))
                scene.add(svgDrw.Text((cxText, cyText), str(nj), text_anchor="middle", size=sizeTextTicks, transform="rotate(90,%s,%s)" % (cxText, cyText)))
                scene.add(svgDrw.Line((offset_matrix_x, cyy), (width-sizeCase, cyy), width=sizeCase*0.1))

    if nx < 300 and ny < 300:
        for i in range(nbLinesX):
            cxLine = i*sizeCase
            scene.add(svgDrw.Line((cxLine+offset_matrix_x, height-offset_matrix_y), (cxLine+offset_matrix_x, height-((nbLinesY-1)*sizeCase+offset_matrix_y)), width=sizeCase*0.01))
        for j in range(nbLinesY):
            cyLine = j*sizeCase
            scene.add(svgDrw.Line((offset_matrix_x, height-(cyLine+offset_matrix_y)), (offset_matrix_x+(nbLinesX-1)*sizeCase, height-(cyLine+offset_matrix_y)), width=sizeCase*0.01))

        chromosome1 = {}
        chromosome2 = {}

        # create chromosomes
        for (i1, s1) in enumerate(genesStrandsC1):
            cx = i1*sizeCase
            symbol = str(symbolsInGenes[0][i1]) if symbolsInGenes is not None else None
            chromosome1[i1] = svgDrw.Gene((cx+offset_matrix_x, height-offset_genes_y), (cx+sizeCase+offset_matrix_x, height-offset_genes_y), strand=s1, width=sizeCase*0.7, stroke_width=0.05*sizeCase, SVGclass=None, text=symbol)
        for (i2, s2) in enumerate(genesStrandsC2):
            cy = i2*sizeCase
            symbol = str(symbolsInGenes[1][i2]) if symbolsInGenes is not None else None
            chromosome2[i2] = svgDrw.Gene((offset_genes_x, height-(cy+offset_matrix_y)), (offset_genes_x, height-(cy+sizeCase+offset_matrix_y)), strand=s2, width=sizeCase*0.7, stroke_width=0.05*sizeCase, SVGclass=None, text=symbol)

        # give a color to each gene using homology relationships
        for (tbs1, tbs2) in homologyGroupsInWindow:
            # Choose a level different from the direct neighbours
            neighboursLevelsOnBothsChrs = reduce(lambda x, y: x | y, [neighboursLevels(chromosome1, i1) for tb1 in tbs1 for i1 in tb1]) |\
                                          reduce(lambda x, y: x | y, [neighboursLevels(chromosome2, i2) for tb2 in tbs2 for i2 in tb2])
            color = closeColorsGenerator.getLevel(differentFrom=neighboursLevelsOnBothsChrs)
            for tb1 in tbs1:
                for i1 in tb1:
                    chromosome1[i1].SVGclass = "HomologGroup%s" % color
            for tb2 in tbs2:
                for i2 in tb2:
                    chromosome2[i2].SVGclass = "HomologGroup%s" % color

        # give grey levels to genes that have no homology in the window

        def giveGreyLevelsTo(chromosome, tbWithNoHomologyInWindow):
            for tb in tbWithNoHomologyInWindow:
                # Choose a level different from the direct neighbours
                nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosome, i) for i in tb])
                grey = farGraysGenerator.getLevel(differentFrom=nLevels)
                for i in tb:
                    chromosome[i].SVGclass = "NoHomologyInWindow%s" % grey
            return chromosome

        chromosome1 = giveGreyLevelsTo(chromosome1, tbWithNoHomologyInWindowC1)
        chromosome2 = giveGreyLevelsTo(chromosome2, tbWithNoHomologyInWindowC2)

        for i1 in genesRemovedDuringFilteringC1:
            chromosome1[i1].SVGclass = "SpeciesSpecificGenes"
        for i2 in genesRemovedDuringFilteringC2:
            chromosome2[i2].SVGclass = "SpeciesSpecificGenes"

        for i1 in chromosome1:
            scene.add(chromosome1[i1])
        for i2 in chromosome2:
            scene.add(chromosome2[i2])

        # fill homologies with +1,-1 or ? or 0
        nonZeroValues = []
        for i1 in hpSigns:
            for i2 in hpSigns[i1]:
                nonZeroValues.append((i1, i2))
                #s = hpSigns[i1][i2][0][1]
                s = hpSigns[i1][i2]
                cx = i1*sizeCase + float(sizeCase)/2
                cy = i2*sizeCase + float(sizeCase)/2
                assert s == +1 or s == -1 or s is None, "s=%s" % s
                assocValue = (("+" if s == +1 else "-") if s is not None else '?')
                scene.add(svgDrw.Text((cx+offset_matrix_x, height-(cy+sizeText*0.16+offset_matrix_y)), assocValue, text_anchor="middle", size=sizeText))
        if nx < 20 and ny < 20:
            for (i1, i2) in itertools.product(range(nx), range(ny)):
                if (i1, i2) not in nonZeroValues:
                    cx = i1*sizeCase + float(sizeCase)/2
                    cy = i2*sizeCase + float(sizeCase)/2
                    scene.add(svgDrw.Text((cx+offset_matrix_x, height-(cy+sizeText*0.16+offset_matrix_y)), "0", text_anchor="middle", size=sizeText, fill=(200, 200, 200), stroke=None))

    else:
        # represent homologies with a black rectangle
        for i1 in hpSigns:
            for i2 in hpSigns[i1]:
                if (i1, i2) not in [dot for diag in diagsIndices for dot in diag]:
                    cx_s = i1*sizeCase
                    cy_s = i2*sizeCase
                    scene.add(svgDrw.Rectangle((cx_s+offset_matrix_x, height-(cy_s+sizeCase+offset_matrix_y)), sizeCase, sizeCase, fill=(0, 0, 0), fill_opacity=0.90))
        print >> sys.stderr, "Warning : some supplementary informations are not displayed because one of the two dimension of the window is > 300"

    if outputFileName is not None:
        scene.write_svg(filename=str(outputFileName))
    #scene.display()
    return scene.strarray()


def test(outFileName):
    scenario = arguments["scenario"]
    if scenario==1 :
        nx=11
        ny=7
        homologySameStrand=[(1,1), (2,2), (3,3), (4,4), (5,5), (8,3), (9,4)]
        homologyOppositeStrand=[]
        diagsIndices=[[(1,1), (2,2), (3,3), (4,4), (5,5)],[(8,3), (9,4)]]
    elif scenario ==2:
        nx=9
        ny=9
        homologySameStrand=[(1,1), (2,2), (3,4), (3,5), (4,3), (5,3), (6,6), (7,7)] # not necessary to mark all for hps
        homologyOppositeStrand=[(3,3), (4,4), (4,5), (5,4), (5,5)]
        diagsIndices=[homologySameStrand + homologyOppositeStrand]
    elif scenario == 3:
        nx=6
        ny=6
        homologySameStrand=[(1,1), (2,2)]
        homologyOppositeStrand=[(3,3), (4,4)]
        diagsIndices=[homologySameStrand, [(3,3)], [(4,4)]]
    elif scenario == 4:
        nx=10
        ny=7
        homologySameStrand=[(3,3)]
        homologyOppositeStrand=[(1,5), (2,4), (4,3), (5,3), (6,3), (7, 2), (8, 1)]
        diagsIndices=[homologySameStrand + homologyOppositeStrand]
    elif scenario == 5:
        nx=9
        ny=9
        homologySameStrand=[(1,1), (2,2), (3,3), (6,6), (7,7)]
        homologyOppositeStrand=[]
        diagsIndices=[homologySameStrand, homologyOppositeStrand]
    elif scenario == 6:
        nx=13
        ny=13
        homologySameStrand=[(1,6), (2,7), (4,8), (5,9), (6,10), (10,8), (11,9)]
        homologyOppositeStrand=[(7,3),(8,2), (3,8)]
    elif scenario == 7:
        # For the paper on the method to extract diagonals
        nx=8
        ny=9
        homologySameStrand=[(1,2), (2,1), (3,3), (4,4), (4,5), (4,6), (5,4), (5,5), (5,6), (6,7)]
        homologyOppositeStrand=[(1,1),(2,2)]
        diagsIndices=[]
    elif scenario == 8:
        # For the paper on the method to extract diagonals
        nx=6
        ny=6
        homologySameStrand=[(1,1), (2,2), (3,3), (4,4)]
        homologyOppositeStrand=[]
        diagsIndices=[]
    elif scenario == 9:
        # For the paper on the method to extract diagonals
        nx=7
        ny=5
        homologySameStrand=[(2,2)]
        homologyOppositeStrand=[(1,3),(3,2),(4,2),(5,1)]
        diagsIndices=[]
    elif scenario == 10:
        # For the paper on the method to extract diagonals
        nx=5
        ny=5
        homologySameStrand=[]
        homologyOppositeStrand=[(1,3),(2,2),(3,1)]
        diagsIndices = []
    elif scenario == 11:
        # For the paper on the method to extract diagonals (merge and probabilities)
        nx=11
        ny=8
        homologyPlusAssociatedValues=[(1,(5,))]
        homologyMinusAssociatedValues=[(4,(4,)),(5,(3,)),(9,(1,))]
        homologyUnknownAssociatedValues=[(7,(6,)),(8,(2,))]
        hpSigns = collections.defaultdict(dict)
        genesStrandsC1=[]
        genesStrandsC2=[]
        for i1 in range(nx):
            genesStrandsC1.append(random.choice([-1,1]))
        for i2 in range(ny):
            genesStrandsC2.append(random.choice([-1,1]))

        for (indexesAVs,aV) in itertools.izip([homologyPlusAssociatedValues, homologyMinusAssociatedValues, homologyUnknownAssociatedValues],[+1,-1,None]):
            for (i1,i2s) in indexesAVs:
                for i2 in i2s:
                    hpSigns[i1][i2]=aV
                    # Determine possible orientations
                    strand = random.choice([-1,+1])
                    if  aV == +1:
                        genesStrandsC1[i1]=strand
                        genesStrandsC2[i2]=strand
                    elif aV == -1:
                        genesStrandsC1[i1]=strand
                        genesStrandsC2[i2]=-strand
                    elif aV == None:
                        b=random.choice([True, False])
                        genesStrandsC1[i1]=None if b else strand
                        genesStrandsC2[i2]=strand if b else None
                    else:
                        raise ValueError()
        begC1,endC1 = 0,nx-1
        begC2,endC2 = 0,ny-1
        tbWithNoHomologyInWindowC1=[0,2,3,6,10]
        tbWithNoHomologyInWindowC2=[0,7]
        homologyGroupsInWindow=[([1],[5]),([4],[4]),([5],[3]),([9],[1]),([7],[6]),([8],[2])]
        diagsIndices=[[(4,4),(5,3),(8,2),(9,1)]]
        genesRemovedDuringFilteringC1=[]
        genesRemovedDuringFilteringC2=[]
        symbolsInGenes=[[1,1,1,1,1,1,1,2,2,1,1],[1,1,2,1,1,1,2,1]]

        drawHomologyMatrix(((begC1,endC1),(begC2,endC2)), (genesStrandsC1, genesStrandsC2), (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2), (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2), hpSigns, homologyGroupsInWindow, diagsIndices, outputFileName=outFileName, maxWidth=100, maxHeight=100 , symbolsInGenes=symbolsInGenes)

if __name__ == '__main__':
    arguments = myTools.checkArgs([("scenario",int)],[("out:FileName",str,"image.svg")],__doc__)
    test(arguments['out:FileName'])
