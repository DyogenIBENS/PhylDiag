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
import copy

import sys
import collections
import itertools
import random
import math
import utils.myTools as myTools
import utils.mySvgDrawer as svgDrw
import utils.myDiags as myDiags
import utils.myMapping as myMapping
from utils.mySvgDrawer import Point
import utils.myLightGenomes as myLightGenomes


# parse the user input (text) for the chromosome range and asses if this query is consistent with the genome data
# g2gtb is a dictionnary to convert indices from gene coordinates into tb coordinates
def parseChrRange(text, genome, g2gtb=None):
    if len(text.split(":")) == 2:
        chr = text.split(":")[0]
        if chr not in genome.keys():
            print >> sys.stderr, "chr %s not in genome" % chr
            raise ValueError
    else:
        print >> sys.stderr, "range not formated as expected : \"chr:deb-fin\""
        raise ValueError

    range = text.split(":")[-1].split("-")
    if len(range) == 2 and range[0].isdigit and (range[1].isdigit() or range[1] == '~'):
        if g2gtb is not None:
            assert isinstance(g2gtb, dict)
            targetIdxG1 = int(range[0]) - 1
            if targetIdxG1 in g2gtb[chr]:
                range[0] = g2gtb[chr][targetIdxG1]
            else:
                print >> sys.stderr, "Warning, gene %s at %s:%s is not part of a tandem block, thus we take the nearest tandem block" %\
                                     (genome[chr][targetIdxG1].n, chr, targetIdxG1 + 1)
                idxGeneThatIsInTargetTb1 = min(g2gtb[chr].keys(), key=lambda x: abs(x-targetIdxG1))
                print >> sys.stderr, "Warning: abs(idxGeneThatIsInTargetTb1 - targetIdxG1) = %s" % abs(idxGeneThatIsInTargetTb1 - targetIdxG1)
                range[0] = g2gtb[chr][idxGeneThatIsInTargetTb1]
            if range[1] == '~':
                range[1] = max(idxTb for idxTb in g2gtb[chr])
            else:
                targetIdxG2 = int(range[1]) - 1
                if targetIdxG2 in g2gtb[chr]:
                    range[1] = g2gtb[chr][targetIdxG2] + 1
                else:
                    print >> sys.stderr, "Warning, gene %s at %s:%s is not part of a tandem block, thus we take the nearest tandem block" %\
                                         (genome[chr][targetIdxG2].n, chr, targetIdxG2 + 1)
                    idxGeneThatIsInTargetTb2 = min(g2gtb[chr].keys(), key=lambda x: abs(x-targetIdxG2))
                    print >> sys.stderr, "Warning: abs(idxGeneThatIsInTargetTb2 - targetIdxG2) = %s" % abs(idxGeneThatIsInTargetTb2 - targetIdxG2)
                    range[1] = g2gtb[chr][idxGeneThatIsInTargetTb2] + 1
        else:
            range = (int(range[0])-1, int(range[1]) if range[1] != '~' else len(genome[chr]))
        if range[0] < 0 or len(genome[chr]) < range[1] or range[1] <= 0 or range[1] <= range[0]:
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
    HomologyGroup1 = set()
    HomologyGroup2 = set()
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


def genesComputeHomologyInformations(chr1, chr2, chrom1, chrom2, families,
                                     filterType,
                                     minChromLength, tandemGapMax):
    c1_fID = myMapping.labelWithFamID(chrom1, families)
    c2_fID = myMapping.labelWithFamID(chrom2, families)
    # Must be applied on the two genomes
    ((c1_fID_filt, Cf2CfID1, (nCL1, nGL1)),
     (c2_fID_filt, Cf2CfID2, (nCL2, nGL2))) = myDiags.filter2D(c1_fID, c2_fID,
                                                               filterType,
                                                               minChromLength,
                                                               keepOriginal=True)
    (chrom1_tb, Ctb2Cf1, nGTD1) = myMapping.remapRewriteInTb(c1_fID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=None)
    (chrom2_tb, Ctb2Cf2, nGTD2) = myMapping.remapRewriteInTb(c2_fID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=None)
    Ctb2CfID1 = {}
    for c in Ctb2Cf1:
        # see Mapping class addition
        Ctb2CfID1[c] = Ctb2Cf1[c] + Cf2CfID1[c]
    Ctb2CfID2 = {}
    for c in Ctb2Cf2:
        # see Mapping class addition
        Ctb2CfID2[c] = Ctb2Cf2[c] + Cf2CfID2[c]

    #Focus on the chromosome of the window
    chrom1_ = chrom1[chr1]
    chrom2_ = chrom2[chr2]
    c1_fID = c1_fID[chr1]
    c2_fID = c2_fID[chr2]
    c1_fID_filt = c1_fID_filt[chr1]
    c2_fID_filt = c2_fID_filt[chr2]
    Cf2CfID1 = Cf2CfID1[chr1]
    Cf2CfID2 = Cf2CfID2[chr2]
    chrom1_tb = chrom1_tb[chr1]
    chrom2_tb = chrom2_tb[chr2]
    Ctb2Cf1 = Ctb2Cf1[chr1]
    Ctb2Cf2 = Ctb2Cf2[chr2]
    Ctb2CfID1 = Ctb2CfID1[chr1]
    Ctb2CfID2 = Ctb2CfID2[chr2]

    ###
    # Build genesRemovedDuringFilteringCX = [..., i, ...] the list of genes that
    # have been removed during the filtering process
    ###
    genesRemovedDuringFilteringC1 = [i1 for (i1, (anc, _)) in enumerate(c1_fID) if anc is None]
    genesRemovedDuringFilteringC2 = [i2 for (i2, (anc, _)) in enumerate(c2_fID) if anc is None]

    (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow) =\
        TbComputeHomologyInformations(chrom1_tb, chrom2_tb)

    genesHomologiesHpSign = collections.defaultdict(lambda: collections.defaultdict(int))
    for i1_tb in TbHpSign:
        for i2_tb in TbHpSign[i1_tb]:
            for (i1, i2) in itertools.product([ii1 for ii1 in Ctb2CfID1[i1_tb]],
                                              [ii2 for ii2 in Ctb2CfID2[i2_tb]]):
                s1 = chrom1_[i1][1]
                s2 = chrom2_[i2][1]
                genesHomologiesHpSign[i1][i2] = s1*s2

    ###
    # Build genesNoHomologiesInWindow1 = [..., [i5,i6,i7], ...] list of tbs with [i5,i6,i7] a tb of three genes whose indices are i5,i6 and i7
    ###
    genesNoHomologiesInWindowC1 = []
    genesNoHomologiesInWindowC2 = []
    for i1_tb in TbNoHomologiesInWindowC1:
        genesNoHomologiesInWindowC1.append([i1 for i1 in Ctb2CfID1[i1_tb]])
    for i2_tb in TbNoHomologiesInWindowC2:
        genesNoHomologiesInWindowC2.append([i2 for i2 in Ctb2CfID2[i2_tb]])

    ###
    # Build genesHomologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tb_CX_X : [..., i, ...] list of indices of genes in the tb
    ###
    genesHomologyGroupsInWindow = []
    for (tbs1, tbs2) in TbHomologyGroupsInWindow:
        genesHomologyGroupsInWindow.append(([], []))
        for i1_tb_group in tbs1:
            genesHomologyGroupsInWindow[-1][0].append([i1 for i1 in Ctb2CfID1[i1_tb_group]])
        for i2_tb_group in tbs2:
            genesHomologyGroupsInWindow[-1][1].append([i2 for i2 in Ctb2CfID2[i2_tb_group]])

    return ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
            genesHomologiesHpSign,
            (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
            genesHomologyGroupsInWindow)


# Generator of levels for colors or gray indices within a palette:
# farIdxs may be an int. The more this int is high, the more neighbour color will be different
class levelIdxGenerator():
    def __init__(self, farIdxs=None, grays=False):
        if grays is not None:
            # see the css joint file "HomologyGroup" ranges from 0 to 44 (included)
            self.firstLevelIdx = 3
            self.lastLevelIdx = 44
        else:
            # see the css joint file "NoHomologyInWindow" ranges from 0 to 14 (included)
            self.firstLevelIdx = 3
            self.lastLevelIdx = 14

        if farIdxs is None:
            self.availableLevels = range(self.firstLevelIdx, self.lastLevelIdx+1)
        else:
            assert farIdxs < self.lastLevelIdx
            self.availableLevels = []
            tmp = range(0, farIdxs)
            random.shuffle(tmp)
            for i in tmp:
                self.availableLevels += range(self.firstLevelIdx + i, self.lastLevelIdx + 1, farIdxs)
            print >> sys.stderr, self.availableLevels
        self.currIdx = 0
        assert len(self.availableLevels) == (self.lastLevelIdx+1) - self.firstLevelIdx

    def getLevel(self, differentFrom=set()):
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
    levelL = chromosome[i-1].SVGclass if 0 <= i-1 <= len(chromosome)-1 else None
    # level of the 'R'ight neighbour
    levelR = chromosome[i+1].SVGclass if 0 <= i+1 <= len(chromosome)-1 else None

    def convertSVGclassIntoInt(strSVGclass):
        if strSVGclass is not None:
            try:
                return int(strSVGclass[-2:])
            except:
                if strSVGclass == 'SpeciesSpecificGenes':
                    return None
                else:
                    return int(strSVGclass[-1])

    levelL = convertSVGclassIntoInt(levelL)
    levelR = convertSVGclassIntoInt(levelR)
    neighboursLevels = {levelL, levelR}
    neighboursLevels = set([l for l in neighboursLevels if l is not None])
    return neighboursLevels


def prepareChromosome(genesStrands, homologousTbs=None, tbsWithNoHomolog=None, genesRemovedDuringFiltering=None,
                      noHomologGraysGenerator=None, homologsColorsGenerator=None, symbolsInGenes=None, lengthGene=1):

    chromosomeItems = []

    # homologousTbs = [homolog1[tb1=[gene1Idx, ...], tb2=[geneAIdx, ...]],
    #                  homolog2[tb1'=[gene1'Idx, ...], tb2'=[geneA'Idx, ...]],
    #                 ... ]
    if not noHomologGraysGenerator:
        noHomologGraysGenerator = levelIdxGenerator(farIdxs=None, grays=True)
    if not homologsColorsGenerator:
        homologsColorsGenerator = levelIdxGenerator(farIdxs=None)

    def giveGreyLevelsTo(chromosome, tbsWithNoHomologyInWindow=None):
        if tbsWithNoHomologyInWindow:
            for tb in tbsWithNoHomologyInWindow:
                # Choose a level different from the direct neighbours
                nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosome, i) for i in tb])
                grey = noHomologGraysGenerator.getLevel(differentFrom=nLevels)
                for i in tb:
                    chromosome[i].SVGclass = "NoHomologyInWindow%s" % grey
        return chromosome

    # create chromosomes
    for (i, s) in enumerate(genesStrands):
        cx = i * lengthGene
        symbol = symbolsInGenes[i] if symbolsInGenes is not None else None
        chromosomeItems.append(svgDrw.Gene(Point(cx, 0),
                                    Point(cx + lengthGene, 0),
                                    strand=s, width=lengthGene*0.7, stroke_width=0.05*lengthGene, SVGclass=None,
                                    text=symbol))

    # give grey levels to genes that have no homology in the window
    chromosomeItems = giveGreyLevelsTo(chromosomeItems, tbsWithNoHomolog)

    if genesRemovedDuringFiltering:
        for i in genesRemovedDuringFiltering:
            chromosomeItems[i].SVGclass = "SpeciesSpecificGenes"

    if homologousTbs:
        for tbs in homologousTbs:
            nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosomeItems, i) for tb in tbs for i in tb])
            color = homologsColorsGenerator.getLevel(differentFrom=nLevels)
            for tb in tbs:
                for i in tb:
                    chromosomeItems[i].SVGclass = "HomologGroup%s" % color
    return chromosomeItems

# TODO, draw in tbs option
def drawChromFromLightGenome(genome, chr, families=None, tandemGapMax=None, lengthGene=1, homologsColorsGenerator=None):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert families is None or isinstance(families, myLightGenomes.Families)

    # FIXME
    newGenome = myLightGenomes.LightGenome()
    newGenome[chr] = genome[chr]
    genome = newGenome

    genesStrands = [s for (_, s) in genome[chr]]

    if families:
        genome_fam = myMapping.labelWithFamNames(genome, families)
        (genome_fam_filt, Cf2Cfam, (nbChrLoss, nbGeneLoss)) = \
            myMapping.remapFilterGeneContent(genome_fam, removedNames={None}, mOld=None)
    else:
        Cf2Cfam = {}
        for c in genome.keys():
            Cf2Cfam[c] = myMapping.Mapping([[i] for i,_ in enumerate(genome[c])])
        genome_fam = genome
    print >> sys.stderr, genome_fam

    if tandemGapMax:
        (genome_tb, Ctb2Cf, nGTD) = myMapping.remapRewriteInTb(genome_fam_filt, tandemGapMax=tandemGapMax, mOld=None)
    else:
        Ctb2Cf = {}
        for c in genome.keys():
            Ctb2Cf[c] = myMapping.Mapping([[i] for i, _ in enumerate(genome_fam[c])])
        genome_tb = genome_fam

    Ctb2Cfam = {}
    for c in Ctb2Cf:
        # see Mapping class addition
        Ctb2Cfam[c] = Ctb2Cf[c] + Cf2Cfam[c]

    ###
    # Build genesRemovedDuringFilteringCX = [..., i, ...] the list of genes that
    # have been removed during the filtering process
    ###
    genesRemovedDuringFiltering = [i1 for (i1, (anc, _)) in enumerate(genome_fam[chr]) if anc is None]

    homologousTbs = []
    genome_tb.computeDictG2Ps()
    setOfPositionsOfTb = {}
    for fam in genome_fam.getGeneNames(checkNoDuplicates=False):
        famPositions = genome_tb.getPositions(fam, default=None)
        if famPositions:
            setOfPositionsOfTb[fam] = famPositions

    for fam, famPositions in setOfPositionsOfTb.iteritems():
        homologousTbs.append([])
        for tbPos in famPositions:
            homologousTbs[-1].append([])
            for i in Ctb2Cfam[tbPos.c][tbPos.idx]:
                homologousTbs[-1][-1].append(i)

    ###
    # Build symbolsInGenes : [ 4,5,1,1,6,2, ...] number of genes in each TB of C
    ###
    symbolsInGenes = [g.n for (i_tb, g) in enumerate(genome_fam[chr])]
    # add a line between genes
    chromosomeItems = []
    chromosomeItems.append(svgDrw.Line(Point(0,0), Point(lengthGene * len(genome[chr]), 0)))
    chromosomeItems.extend(prepareChromosome(genesStrands, homologousTbs=homologousTbs,
                                   tbsWithNoHomolog=None,
                                   genesRemovedDuringFiltering=genesRemovedDuringFiltering,
                                   noHomologGraysGenerator=None,
                                   homologsColorsGenerator=homologsColorsGenerator,
                                   symbolsInGenes=symbolsInGenes,
                                   lengthGene=lengthGene))
    return chromosomeItems

def drawChromosomes(genesStrandsC1, tbWithNoHomologyInWindowC1, genesRemovedDuringFilteringC1,
                    genesStrandsC2, tbWithNoHomologyInWindowC2, genesRemovedDuringFilteringC2,
                    homologyGroupsInWindow, closeColorsGenerator,
                    symbolsInGenes, sizeCase, height):

    chromosome1 = prepareChromosome(genesStrandsC1, None,  tbWithNoHomologyInWindowC1, genesRemovedDuringFilteringC1,
                                    symbolsInGenes=symbolsInGenes[0] if symbolsInGenes else None,
                                    lengthGene=sizeCase)

    chromosome2 = prepareChromosome(genesStrandsC2, None, tbWithNoHomologyInWindowC2, genesRemovedDuringFilteringC2,
                                    symbolsInGenes=symbolsInGenes[1] if symbolsInGenes else None,
                                    lengthGene=sizeCase)

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

    # insert at the beginning to draw first
    chromosome1.insert(0, svgDrw.Line(Point(0, 0), Point(len(genesStrandsC1) * sizeCase, 0)))
    chromosome2.insert(0, svgDrw.Line(Point(0, 0), Point(len(genesStrandsC2) * sizeCase, 0)))
    return (chromosome1, chromosome2)

def drawMatrix(nx, ny, (begC1, endC1), (begC2, endC2), hpSigns, diagsIndices, sizeCell, width, height,
               diagColorGenerator=None, scaleFactorRectangles=1.0):
        print >> sys.stderr, scaleFactorRectangles
        assert isinstance(scaleFactorRectangles, float)
        sizeText = float(sizeCell*0.9)
        listOfMatrixItems = []
        if not diagColorGenerator:
            diagColorGenerator = levelIdxGenerator(farIdxs=True)

        nbLinesX = nx+1 #Nb of vertical lines (x varies) in the matrix
        nbLinesY = ny+1 #Nb of horizontal lines (y varies) in the matrix

        # draw Diagonals first because they are on the background
        print >> sys.stderr, "Nb of diagonals showed = ", len(diagsIndices)

        # sort diagonals to give colors according to localisations of diagonals
        diagsIndices.sort(key=lambda x: x[0])
        for diag in diagsIndices:
            # choose a color different from the neighbours
            color = diagColorGenerator.getLevel()
            for (i, j) in diag:
                cx_s = i*sizeCell
                cy_s = j*sizeCell
                if nx >= 300 or ny >= 300:
                    listOfMatrixItems.append(svgDrw.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                               sizeCell, sizeCell, fill_opacity=0.90, svgClass="HomologGroup%s" % color,
                                               bigger=scaleFactorRectangles))
                else:
                    listOfMatrixItems.append(svgDrw.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                               sizeCell, sizeCell, fill_opacity=0.90, svgClass="HomologGroup%s" % color))
                #sizeCell(scgDrw.Rectangle((cx_s,height-(cy_s+sizeCell)), sizeCell, sizeCell, fill_opacity=0.90, svgClass = "chromgrp%s" % color ))
                #listOfMatrixItems.append(svgDrw.Rectangle((cx_s,height-(cy_s+sizeCell)), sizeCell, sizeCell, fill_opacity=fill_opacity=0.90, color=(color,color,color))
            # draw rectangles around diagonals
            min_i = min(diag, key=lambda x: x[0])[0]
            max_i = max(diag, key=lambda x: x[0])[0]
            min_j = min(diag, key=lambda x: x[1])[1]
            max_j = max(diag, key=lambda x: x[1])[1]
            cx_s = min_i*sizeCell
            cy_s = max_j*sizeCell

            scaleFactorBoundingBoxDiags = 1.5 * scaleFactorRectangles if nx >= 300 or ny >= 300 else 1.0
            listOfMatrixItems.append(svgDrw.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                       (max_j-min_j)*sizeCell + sizeCell, (max_i-min_i)*sizeCell + sizeCell,
                                       stroke='black', fill='none', strokeWidth=0.2*sizeCell*scaleFactorBoundingBoxDiags))

        # tick lines
        widthTicks = min(float(width)/1000, float(height)/1000)
        sizeTextTicks = widthTicks*10
        for (i, ni) in enumerate(range(begC1, endC1)):
            cx = i*sizeCell
            if ni % 10 == 0:
                # ticks
                listOfMatrixItems.append(svgDrw.Line(Point(sizeCell/2+cx, height + sizeCell/2),
                                      Point(sizeCell/2+cx, height), width=widthTicks))
            if ni % 50 == 0:
                cyText = height - max(sizeCell/2, sizeTextTicks/2)
                cxx = cx + sizeCell/2
                if nx > 750 or ny > 750:
                    if ni % 100 == 0:
                        # TODO FIXME, too high for long chromosomes
                        listOfMatrixItems.append(svgDrw.Text(Point(cxx, cyText + sizeCell), str(ni), text_anchor="middle", size=sizeTextTicks))
                        #listOfMatrixItems.append(svgDrw.Text(Point(cxx, cyText), str(ni), text_anchor="middle", size=sizeTextTicks))
                        listOfMatrixItems.append(svgDrw.Line(Point(cxx, ny * sizeCell),
                                              Point(cxx, sizeCell), width=sizeCell*0.1))
                else:
                    listOfMatrixItems.append(svgDrw.Text(Point(cxx, cyText + sizeCell), str(ni), text_anchor="middle", size=sizeTextTicks))
                    listOfMatrixItems.append(svgDrw.Line(Point(cxx, height),
                                          Point(cxx, height - ny * sizeCell), width=sizeCell*0.1))

        for (j, nj) in enumerate(range(begC2, endC2)):
            cy = j*sizeCell
            if nj % 10 == 0:
                listOfMatrixItems.append(svgDrw.Line(Point(sizeCell/2 - sizeCell, height - (sizeCell/2+cy)),
                                      Point(0, (height - (sizeCell/2+cy))), width=widthTicks))
            if nj % 50 == 0:
                cxText = - max(sizeCell/2, sizeTextTicks/2)
                cyy = height - (sizeCell/2 + cy)
                if nx > 750 or ny > 750:
                    if nj % 100 == 0:
                        # TODO FIXME, not visible for big chromosomes
                        #listOfMatrixItems.append(svgDrw.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCell,cxx+sizeCell,cyy)))
                        listOfMatrixItems.append(svgDrw.Text(Point(0, 0),  str(nj),
                                          text_anchor="middle", size=sizeTextTicks, transform="translate(%s,%s) rotate(-90)" % (cxText, cyy - 6*sizeCell)))
                        # listOfMatrixItems.append(svgDrw.Text(Point(cxText, cyText), str(nj),
                        #                       text_anchor="middle", size=sizeTextTicks, transform="rotate(90,%s,%s)" % (cxText, cyText)))
                        listOfMatrixItems.append(svgDrw.Line(Point(0, cyy),
                                              Point(width-sizeCell, cyy), width=sizeCell*0.1))
                else:
                    #listOfMatrixItems.append(svgDrw.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCell,cxx+sizeCell,cyy)))
                    # FIXME: why - 6*sizeCell ?
                    listOfMatrixItems.append(svgDrw.Text(Point(0, 0),  str(nj),
                                          text_anchor="middle", size=sizeTextTicks, transform="translate(%s,%s) rotate(90)" % (cxText, cyy - 6*sizeCell)))
                    listOfMatrixItems.append(svgDrw.Line(Point(0, cyy),
                                          Point(nx * sizeCell, cyy), width=sizeCell*0.1))

        if nx < 300 and ny < 300:
            for i in range(nbLinesX):
                cxLine = i*sizeCell
                listOfMatrixItems.append(svgDrw.Line(Point(cxLine, height),
                                      Point(cxLine, height-((nbLinesY-1)*sizeCell)), width=sizeCell*0.01))
            for j in range(nbLinesY):
                cyLine = j*sizeCell
                listOfMatrixItems.append(svgDrw.Line(Point(0, height-(cyLine)),
                                      Point((nbLinesX-1)*sizeCell, height-(cyLine)), width=sizeCell*0.01))

            # fill homologies with +1, -1 or ? or 0
            nonZeroValues = []
            for i1 in hpSigns:
                for i2 in hpSigns[i1]:
                    nonZeroValues.append((i1, i2))
                    #s = hpSigns[i1][i2][0][1]
                    s = hpSigns[i1][i2]
                    cx = i1*sizeCell + float(sizeCell)/2
                    cy = i2*sizeCell + float(sizeCell)/2
                    assert s == +1 or s == -1 or s is None, "s=%s" % s
                    assocValue = (("+" if s == +1 else "-") if s is not None else '?')
                    listOfMatrixItems.append(svgDrw.Text(Point(cx, height-(cy+sizeText*0.16)),
                                          assocValue, text_anchor="middle", size=sizeText))
            if nx < 20 and ny < 20:
                for (i1, i2) in itertools.product(range(nx), range(ny)):
                    if (i1, i2) not in nonZeroValues:
                        cx = i1*sizeCell + float(sizeCell)/2
                        cy = i2*sizeCell + float(sizeCell)/2
                        listOfMatrixItems.append(svgDrw.Text(Point(cx, height-(cy+sizeText*0.16)),
                                              "0", text_anchor="middle", size=sizeText, fill=(200, 200, 200), stroke=None))
        else:
            # represent homologies with a black rectangle
            for i1 in hpSigns:
                for i2 in hpSigns[i1]:
                    if (i1, i2) not in [dot for diag in diagsIndices for dot in diag]:
                        cx_s = i1*sizeCell
                        cy_s = i2*sizeCell
                        listOfMatrixItems.append(svgDrw.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                                   sizeCell, sizeCell, fill=(0, 0, 0), fill_opacity=0.90, bigger=scaleFactorRectangles))
            print >> sys.stderr, "Warning : some supplementary informations are not displayed because one of the two dimension of the window is > 300"

        return listOfMatrixItems


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
                       outputFileName=None, maxWidth=100, maxHeight=100, symbolsInGenes=None, scaleFactorRectangles=1):
    # example
    # genesStrandsC1 = [1,-1, None, ...]
    # tbWithNoHomologyInWindowC1=[0,2,3,6,10]
    # tbWithNoHomologyInWindowC2=[0,7]
    # homologyGroupsInWindow=[([1],[5]),([4],[4]),([5],[3]),([9],[1]),([7],[6]),([8],[2])]
    # diagsIndices=[[(4,4),(5,3),(8,2),(9,1)]]
    # genesRemovedDuringFilteringC1=[]
    # genesRemovedDuringFilteringC2=[]
    # symbolsInGenes=[[1,1,1,1,1,1,1,2,2,1,1],[1,1,2,1,1,1,2,1]]
    # hpSigns a collections.defaultdict() with hpSigns[i1][i2]=aV
    assert isinstance(genesStrandsC1, list) and isinstance(genesStrandsC2, list)
    assert isinstance(genesRemovedDuringFilteringC1, list) and isinstance(genesRemovedDuringFilteringC2, list)
    assert isinstance(tbWithNoHomologyInWindowC1, list) and isinstance(tbWithNoHomologyInWindowC2, list)
    assert isinstance(homologyGroupsInWindow, list)
    assert isinstance(hpSigns, dict)
    assert isinstance(diagsIndices, list)
    # For the print on screen
    begC1 = begC1 + 1
    begC2 = begC2 + 1
    endC1 = endC1 + 1
    endC2 = endC2 + 1
    # nx = number of genes on the first genome
    # ny = number of genes on the second genome
    nx = len(genesStrandsC1)
    ny = len(genesStrandsC2)

    # the size of the components of the matrix is chosen using the smallest and more restricting dimension
    # (contains more genes comparing to its size)
    sizeCase = float(min(float(maxWidth) / (nx + 3), float(maxHeight) / (ny + 3)))  # +3 for margins

    width = float(nx+3) * sizeCase
    height = float(ny+3) * sizeCase
    print >> sys.stderr, "width=", width
    print >> sys.stderr, "height=", height

    scene = svgDrw.Scene(name='homology_matrix', width=width, height=height)

    closeColorsGenerator = levelIdxGenerator()

    # draw lines of chromosomes
    # offset_genes : corresponds to the chromosome lines positions passing through the middle of the genes
    offset_genes_x = sizeCase + sizeCase/2
    offset_genes_y = sizeCase + sizeCase/2
    # scene.add(svgDrw.Line(Point(offset_genes_x, height - offset_genes_y),
    #                       Point(width-sizeCase, height - offset_genes_y), width=0.1*sizeCase))
    # scene.add(svgDrw.Line(Point(offset_genes_x, sizeCase),
    #                       Point(offset_genes_x, height - (offset_genes_y)), width=0.1*sizeCase))
    offset_matrix_x = 2 * sizeCase
    offset_matrix_y = 2 * sizeCase

    listOfMatrixItems = drawMatrix(nx, ny, (begC1, endC1), (begC2, endC2), hpSigns, diagsIndices, sizeCase, width, height,
                         diagColorGenerator=None, scaleFactorRectangles=scaleFactorRectangles)
    listOfMatrixItems = svgDrw.tanslateItems(listOfMatrixItems, 2 * sizeCase, - 2 * sizeCase)

    listOfItems = []
    if nx < 300 and ny < 300:
        (chromosome1, chromosome2) = drawChromosomes(genesStrandsC1, tbWithNoHomologyInWindowC1, genesRemovedDuringFilteringC1,
                                                     genesStrandsC2, tbWithNoHomologyInWindowC2, genesRemovedDuringFilteringC2,
                                                     homologyGroupsInWindow, closeColorsGenerator,
                                                     symbolsInGenes, sizeCase,
                                                     # offset_genes_x, offset_genes_y,
                                                     # offset_matrix_x, offset_matrix_y,
                                                     height)

        chromosome1 = svgDrw.tanslateItems(chromosome1, 2 * sizeCase, height)
        for item in chromosome2:
            if isinstance(item, svgDrw.Gene):
                item.start = Point(0, height - item.start.x - 2 * sizeCase)
                item.end = Point(0, height - item.end.x - 2 * sizeCase)
            elif isinstance(item, svgDrw.Line):
                item.start = Point(0, height - item.start.x - 2 * sizeCase)
                item.end = Point(0, height - item.end.x - 2 * sizeCase)
            else:
                raise TypeError
        listOfItems = chromosome1 + chromosome2

    listOfItems += listOfMatrixItems
    listOfItems = svgDrw.tanslateItems(listOfItems, sizeCase, -sizeCase)

    for item in listOfItems:
        scene.add(item)
    if outputFileName is not None:
        scene.write_svg(filename=str(outputFileName))
    #scene.display()
    return scene.strarray()

def test(arguments):
    scenario = arguments["scenario"]
    outFileName = arguments['out:fileName']
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
        begC1, endC1 = 0, nx-1
        begC2, endC2 = 0, ny-1
        tbWithNoHomologyInWindowC1=[0,2,3,6,10]
        tbWithNoHomologyInWindowC2=[0,7]
        homologyGroupsInWindow=[([1],[5]),([4],[4]),([5],[3]),([9],[1]),([7],[6]),([8],[2])]
        diagsIndices=[[(4,4),(5,3),(8,2),(9,1)]]
        genesRemovedDuringFilteringC1=[]
        genesRemovedDuringFilteringC2=[]
        symbolsInGenes=[[1,1,1,1,1,1,1,2,2,1,1],[1,1,2,1,1,1,2,1]]

        drawHomologyMatrix(((begC1, endC1),(begC2, endC2)), (genesStrandsC1, genesStrandsC2),
                           (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
                           (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2),
                           hpSigns, homologyGroupsInWindow, diagsIndices,
                           outputFileName=outFileName, maxWidth=100, maxHeight=100 ,
                           symbolsInGenes=symbolsInGenes)

    elif scenario == 12:
        genome = myLightGenomes.LightGenome()
        OG = myLightGenomes.OGene
        genome['0'] = [OG('a', +1), OG('b', +1), OG('c', -1), OG('d', +1), OG('e', -1), OG('f', +1), OG('g', +1), OG('h', +1)]
        lengthGene = 1
        width = len(genome['0']) + 2 * lengthGene
        height = 3 * lengthGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        families = myLightGenomes.Families()
        F = myLightGenomes.Family
        families.addFamily(F('A', ['a', 'b']))
        families.addFamily(F('B', ['c']))
        families.addFamily(F('C', ['d', 'e']))
        families.addFamily(F('D', ['f', 'g', 'h']))
        chromosomeItems = drawChromFromLightGenome(genome, '0', families, lengthGene=lengthGene)
        for item in chromosomeItems:
            # if isinstance(item, svgDrw.Gene):
            # g = item
            item.start = Point(item.start.x, item.start.y + lengthGene)
            item.end = Point(item.end.x, item.end.y + lengthGene)
            scene.add(item)
        scene.write_svg(filename=outFileName)
        print >> sys.stderr, "finished!"
    elif scenario == 13:
        families = myLightGenomes.Families()
        F = myLightGenomes.Family
        familiesDict = collections.defaultdict(list)
        allFamilies = 'ABCDEFGHIJK'
        for familyName in allFamilies:
            familiesDict[familyName] = []
        genome = myLightGenomes.LightGenome()
        OG = myLightGenomes.OGene
        print familiesDict
        for geneName in 'abcdefghijklmnopqrstuvw':
            genome['0'].append(OG(geneName, random.choice([+1, -1])))
            famName = random.choice(allFamilies)
            familiesDict[famName].append(geneName)
        for famName, dn in familiesDict.iteritems():
            families.addFamily(F(famName, dn))

        lengthGene = 1
        width = len(genome['0']) + 2 * lengthGene
        height = 3 * lengthGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        chromosomeItems = drawChromFromLightGenome(genome, '0', families, lengthGene=lengthGene)
        for item in chromosomeItems:
            # if isinstance(item, svgDrw.Gene):
            # g = item
            item.start = Point(item.start.x, item.start.y + lengthGene)
            item.end = Point(item.end.x, item.end.y + lengthGene)
            scene.add(item)
        scene.write_svg(filename=outFileName)
        print >> sys.stderr, "finished!"
    elif scenario == 14:
        os.chdir('/home/jlucas/Libs/MagSimus')
        genome1 = myLightGenomes.LightGenome('data/genesST.Mus.musculus.list.bz2')
        genome2 = myLightGenomes.LightGenome('data/genesST.Gallus.gallus.list.bz2')
        families = myLightGenomes.Families('data/ancGenes.Amniota.list.bz2')


        genome_Mouse = myLightGenomes.LightGenome()
        genome_Chicken = myLightGenomes.LightGenome()
        genome_Mouse['5'] = genome1['5'][597-1:627-1]
        genome_Chicken['4'] = genome2['4'][582-1:611-1]

        newFamilies = myLightGenomes.Families()
        old2newFamily = {}
        allFamilies = []
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789':
            allFamilies.append(c)

        familiesInBothGenomes = genome_Mouse.getOwnedFamilyNames(families) & genome_Chicken.getOwnedFamilyNames(families)
        for gene in genome_Mouse['5'] + genome_Chicken['4']:
            oldFamily = families.getFamilyByName(gene.n, default=None)
            if oldFamily:
                if oldFamily.fn in familiesInBothGenomes:
                    if oldFamily.fn not in old2newFamily:
                        newFamily = allFamilies.pop(0)
                        newFamilies.addFamily(myLightGenomes.Family(newFamily, oldFamily.dns))
                        old2newFamily[oldFamily.fn] = newFamily
        families = newFamilies

        lengthGene = 1
        width = max(len(genome_Mouse['5']), len(genome_Chicken['4'])) + 2 * lengthGene
        height = 4 * lengthGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        chromosomeMouseItems = drawChromFromLightGenome(genome_Mouse, '5', families, lengthGene=lengthGene)
        chromosomeChickenItems = drawChromFromLightGenome(genome_Chicken, '4', families, lengthGene=lengthGene)
        for item in chromosomeMouseItems:
            item.start = Point(item.start.x, item.start.y + lengthGene)
            item.end = Point(item.end.x, item.end.y + lengthGene)
            scene.add(item)
        for item in chromosomeChickenItems:
            item.start = Point(item.start.x, item.start.y + 2*lengthGene)
            item.end = Point(item.end.x, item.end.y + 2*lengthGene)
            scene.add(item)
        scene.write_svg(filename=outFileName)

    elif scenario == 15:
        # presentation of the different events
        OG = myLightGenomes.OGene
        genome = myLightGenomes.LightGenome()
        genome['1'] = [OG('A', +1), OG('B', -1), OG('C', +1)]
        genome['2'] = [OG('D', -1), OG('E', +1), OG('F', -1)]

        sizeGene = 1
        listOfChromosomes = []
        for c in genome:
            listOfChromosomes.append(drawChromFromLightGenome(genome, c, lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        import libs.myEvents as mE

        # Gene events
        # tandem duplication
        newGenome = mE.performInsertNewGene(genome, ('B.b', '1', 2, -1), keepOriginalGenome=True)
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        # gene loss
        newGenome = mE.performGeneLoss(genome, ('1', 1), keepOriginalGenome=True)
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        # de novo gene birth
        newGenome = mE.performInsertNewGene(genome, ('G', '1', 2, -1), keepOriginalGenome=True)
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))

        # Chromosomal rearrangements
        # fission
        newGenome = mE.performFission(genome, ('1', 1), keepOriginalGenome=True)
        assert '0' in newGenome.keys()
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '0', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        # fusion
        newGenome = mE.performFusion(genome, (('1', +1), ('2', +1)), keepOriginalGenome=True)
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        # inversion
        newGenome = mE.performInversion(genome, ('1', 1, 3), keepOriginalGenome=True)
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        # reciprocal translocation
        newGenome = mE.performReciprocalTranslocation(genome, (('1', 1, +1), ('2', 1, +1)), keepOriginalGenome=True)
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))
        listOfChromosomes.append(drawChromFromLightGenome(newGenome, '2', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))

        genome['3'] = [OG('A', +1), OG('B', -1), OG('C', +1), OG('D', +1), OG('E', -1), OG('F', +1), OG('G', +1), OG('H', -1), OG('I', +1)]
        listOfChromosomes.append(drawChromFromLightGenome(genome, '3', lengthGene=sizeGene, homologsColorsGenerator=levelIdxGenerator(farIdxs=5)))

        width = (2 + 8) * sizeGene
        height = (2 + len(listOfChromosomes)) * sizeGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        for i, chromosome in enumerate(listOfChromosomes):
            svgDrw.tanslateItems(chromosome, 0, sizeGene + i * sizeGene)
            for item in chromosome:
                scene.add(item)
        scene.write_svg(filename=outFileName)



if __name__ == '__main__':
    #arguments = myTools.checkArgs([("scenario",int)],[("out:FileName",str,"image.svg")],__doc__)
    import os
    import sys
    print sys.stderr, sys.argv
    os.chdir('/home/jlucas/Libs/MagSimus')
    arguments = {}
    arguments['scenario'] = int(sys.argv[1])
    arguments['out:fileName'] = 'toto.svg'
    test(arguments)
