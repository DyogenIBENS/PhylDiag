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
import collections

import utils.myGenomes as myGenomes
import utils.myTools as myTools
# PhylDiag core algorithm
import utils.myDiags as myDiags
import utils.myProbas as myProbas
import utils.myMapping as myMapping

import drawHomologyMatrixWithSBs

__doc__= """
        Show the homology matrix with coloured synteny blocks (also called diagonals).
        - Each colour represents a synteny block.
        - On the x-axis are the genes of the 1st genome in the desired window
        - On the y-axis are the genes of the 2nd genome in the desired window
        - Each coloured rectangle in the matrix represents a filiation relationship restricted to the ancGene species. Couloured rectangle means the corresponding horizontal and vertical genes come from the same ancestral gene in the ancGene species. Take care that coulour genes are not homology relationship. For instance if the ancestor had two paralogs in its genome, there will be two distinct filiation relationships for genes herited from these genes. Genes are coloured if they are in the same filiation relationship class.
        - '+' indicates filiation relationships that have horizontal gene and vertical gene in the same direction on both chromosomes
          '-' indicates filiation relationships that have horizontal gene and vertical gene in opposite directions on each chromosomes
        """

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
    if len(range) == 2 and range[0].isdigit and (range[1].isdigit() or range[1]=='~'):
        if g2gtb != None:
            range[0] = g2gtb[chr][int(range[0])-1]
            range[1] = g2gtb[chr][int(range[1])-1] if range[1]!='~' else len(genome[chr])
        else:
            range = (int(range[0])-1,int(range[1]) if range[1]!='~' else len(genome[chr]) )
        if range[0] < 0 or range[1] > len(genome[chr]) or  range[1] <= 0  or range[1] <= range[0]:
            print >> sys.stderr,\
                "range %s is incoherent for chr %s. FYI chr %s contains %s elements. Be sure that beginning < end of range" % ([range[0]+1, range[1]], chr, chr, len(genome[chr]))
            raise ValueError
    else:
        raise ValueError
    return (chr, range)

def TbComputeHomologyInformations(chrom1_tb, chrom2_tb):

    locC2 = collections.defaultdict(list)
    for (i2_tb,(ianc,_)) in enumerate(chrom2_tb):
        if ianc != -1:
            locC2[ianc].append( i2_tb )
    ###
    # Build homologiesHpSigns = { i1_tb : {i2_tb : hpSign, ...}...} with hpSign
    # the hpSign (s1*s2) of the homology at the (i1,i2) coordinate
    ###
    TbHpSign = collections.defaultdict(lambda:collections.defaultdict(tuple))
    TbHpSign2 =  collections.defaultdict(int)
    for (i1_tb,(ianc,_)) in enumerate(chrom1_tb):
        if ianc !=-1:
            if ianc in locC2: #TODO remove by advance ianc == -1 from chrom1_tb
                for i2_tb in locC2[ianc]:
                    TbHpSign2[i2_tb] = None # Just to perform a search after
                    s1 = chrom1_tb[i1_tb][1]
                    s2 = chrom2_tb[i2_tb][1]
                    TbHpSign[i1_tb][i2_tb] = myDiags.strandProduct(s1,s2)
    ###
    # Build TBNoHomologiesInWindowX = [..., i_tb, ... ] with i_tb the index of
    # the tb of chromosomeX_tb with no homology in the window
    ###
    TbNoHomologiesInWindowC1 = []
    for (i1_tb,_) in enumerate(chrom1_tb):
        if i1_tb not in TbHpSign:
            TbNoHomologiesInWindowC1.append(i1_tb)
    TbNoHomologiesInWindowC2 = []
    for (i2_tb,_) in enumerate(chrom2_tb):
        if i2_tb not in TbHpSign2:
            TbNoHomologiesInWindowC2.append(i2_tb)
    ###
    # Build homologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tbCX_X : [..., i, ...] list of indices of genes in the tb
    ###
    locC1_tbIndx = collections.defaultdict(list)
    locC2_tbIndx = collections.defaultdict(list)
    TbHomologyGroupsInWindow = []
    listOfHPsCoordinates=[]
    for i1_tb in TbHpSign:
        for i2_tb in TbHpSign[i1_tb]:
            listOfHPsCoordinates.append((i1_tb,i2_tb))
    for (i1_tb,i2_tb) in listOfHPsCoordinates:
        locC2_tbIndx[i1_tb].append(i2_tb)
        locC1_tbIndx[i2_tb].append(i1_tb)
    HomologyGroup1={}
    HomologyGroup2={}
    for (i1_tb,i2_tb) in listOfHPsCoordinates:
        if i1_tb not in HomologyGroup1 or i2_tb not in HomologyGroup2:
            TbHomologyGroupsInWindow.append(([i1_tb],[i2_tb]))
            HomologyGroup1[i1_tb] = None
            for i2_tbb in locC2_tbIndx[i1_tb]:
                TbHomologyGroupsInWindow[-1][1].append(i2_tbb)
                HomologyGroup2[i2_tbb] = None
                for i1_tbb in locC1_tbIndx[i2_tbb]:
                    TbHomologyGroupsInWindow[-1][0].append(i1_tbb)
                    HomologyGroup1[i1_tbb] = None
                del locC1_tbIndx[i2_tbb]
            del locC2_tbIndx[i1_tb]
    del HomologyGroup1
    del HomologyGroup2

    return (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow)

def TbComputeDiagIndices(TbListOfDiags):
    ###
    # Build TbDiagIndices = [..., [...,(i16,j16),...], ...] list of diagonals
    # with diagonals = list of all the points of the diagonal
    ###
    TbDiagIndices = []
    for TbDiag in TbListOfDiags:
        ((c1,tbs1),(c2,tbs2),tbsA,pVal) = TbDiag
        TbDiagIndices.append([])
        for (indx,(_,_,_)) in enumerate(tbsA):
            TbDiagIndices[-1].append((tbs1[indx],tbs2[indx]))
    return TbDiagIndices

def genesComputeHomologyInformations(chrom1, chrom2):
    c1_aID = myMapping.labelWithAncGeneID(chrom1, ancGenes)
    c2_aID = myMapping.labelWithAncGeneID(chrom2, ancGenes)
    # Must be applied on the two genomes
    ((c1_aID_filt, Cf2CaID1, (nCL1, nGL1)),
     (c2_aID_filt, Cf2CaID2, (nCL2, nGL2))) =\
        myDiags.filter2D(c1_aID, c2_aID,
                         filterType,
                         arguments["minChromLength"],
                         keepOriginal=True)
    (chrom1_tb, Ctb2Cf1, nGTD1) =\
        myMapping.remapRewriteInTb(c1_aID_filt,
                                   tandemGapMax=arguments['tandemGapMax'],
                                   mOld=None)
    (chrom2_tb, Ctb2Cf2, nGTD2) =\
        myMapping.remapRewriteInTb(c2_aID_filt,
                                   tandemGapMax=arguments['tandemGapMax'],
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
    chrom1 = chrom1[chr1]
    chrom2 = chrom2[chr2]
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
    genesRemovedDuringFilteringC1 = [i1 for (i1,(anc,_)) in enumerate(c1_aID) if anc == None]
    genesRemovedDuringFilteringC2 = [i2 for (i2,(anc,_)) in enumerate(c2_aID) if anc == None]


    (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow) =\
        TbComputeHomologyInformations(chrom1_tb, chrom2_tb)

    genesHomologiesHpSign = collections.defaultdict(lambda:collections.defaultdict(int))
    for i1_tb in TbHpSign:
        for i2_tb in TbHpSign[i1_tb]:
            # TODO : propose to vizualize the mh with the species specific genes
            for (i1,i2) in itertools.product(
                [ii1 for ii1 in Ctb2CaID1[i1_tb]],
                [ii2 for ii2 in Ctb2CaID2[i2_tb]]):
                s1 = chrom1[i1][1]
                s2 = chrom2[i2][1]
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
    for (tbs1,tbs2) in TbHomologyGroupsInWindow:
        genesHomologyGroupsInWindow.append(([],[]))
        for i1_tb_group in tbs1:
            genesHomologyGroupsInWindow[-1][0].append([i1 for i1 in Ctb2CaID1[i1_tb_group]])
        for i2_tb_group in tbs2:
            genesHomologyGroupsInWindow[-1][1].append([i2 for i2 in Ctb2CaID2[i2_tb_group]])

    return ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
            genesHomologiesHpSign,
            (genesNoHomologiesInWindowC1,genesNoHomologiesInWindowC2),
            genesHomologyGroupsInWindow)

def genesComputeDiagIndices(listOfDiags):
    ###
    # Build diagIndices = [..., [...,(i16,j16),...], ...] mist of diagonals with
    #   diagonals = list of all the points of the diagonal
    ###
    diagIndices = []
    for diag in listOfDiags:
        ((c1,l1), (c2,l2), la, pVal) = diag
        diagIndices.append([])
        for (idxHp, (anc, strand, dist)) in enumerate(la):
            diagIndices[-1].extend(itertools.product(l1[idxHp], l2[idxHp]))
    return diagIndices

arguments = myTools.checkArgs(
    [("genome1", file), ("genome2", file),
     ("ancGenes", file), ("chr1:deb1-fin1", str),
     ("chr2:deb2-fin2", str)],
    [("gapMax", str, 'None'),
     ("tandemGapMax", int, 0),
     ("consistentSwDType", bool, True),
     ("filterType", str, 'InCommonAncestor'),
     ("minChromLength", int, 1),
     ("pThreshold", float, 0.001),
     ("out:SyntenyBlocks", str, "./res/syntenyBlocksDrawer.txt"),
     ("mode:chromosomesRewrittenInTbs", bool, False),
     ('convertGenicToTbCoordinates', bool, False),
     ("distanceMetric", str, 'CD'),
     ('nbHpsRecommendedGap', int, 2),
     ('targetProbaRecommendedGap', float, 0.01),
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

assert arguments["distanceMetric"] == 'DPD' or arguments["distanceMetric"] == 'MD'\
    or arguments["distanceMetric"] == 'CD' or arguments["distanceMetric"] == 'ED'
assert (arguments["convertGenicToTbCoordinates"] and
        arguments["mode:chromosomesRewrittenInTbs"])\
        or not arguments["convertGenicToTbCoordinates"]

# Load genomes
genome1 = myGenomes.Genome(arguments["genome1"])
genome2 = myGenomes.Genome(arguments["genome2"])
# Change genome format
genome1Name = genome1.name
genome2Name = genome2.name
genome1_ = {}
for chrom in genome1.lstGenes:
    genome1_[chrom] = genome1.lstGenes[chrom]
    genome1_[chrom] = [(g.names[0],g.strand) for g in genome1_[chrom]]
genome2_ ={}
for chrom in genome2.lstGenes:
    genome2_[chrom] = genome2.lstGenes[chrom]
    genome2_[chrom] = [(g.names[0],g.strand) for g in genome2_[chrom]]
genome1=genome1_
genome2=genome2_
ancGenes = myGenomes.Genome(arguments["ancGenes"])
modesFilter = list(myDiags.FilterType._keys)
filterType = myDiags.FilterType[modesFilter.index(arguments["filterType"])]
thresholdChr = 50

if not arguments['mode:chromosomesRewrittenInTbs']:
    #chromosomes are shown as a list of genes
    #print >> sys.stderr, "List of (chromosomes, length in genes) of Genome 1, for chr of size > %s " % thresholdChr
    #for (chr1,len1) in [(key1, len(chr1)) for (key1,chr1) in genome1.items() if len(chr1) > thresholdChr]:
    #       print >> sys.stderr, "chr %s has %s genes" % (chr1,len1)
    #print >> sys.stderr, "List of (chromosomes, length in genes) of Genome 2, for chr of size > %s " % thresholdChr
    #for (chr2,len2) in [(key2, len(chr2)) for (key2,chr2) in genome2.items() if len(chr2) > thresholdChr]:
    #       print >> sys.stderr, "chr %s has %s genes" % (chr2,len2)
    (chr1,range1) = parseChrRange(arguments["chr1:deb1-fin1"], genome1)
    (chr2,range2) = parseChrRange(arguments["chr2:deb2-fin2"], genome2)
    chrom1 ={}
    chrom2 ={}
    chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
    chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]

    ###
    # Build Genes Strands
    ###
    genesStrandsC1 = [s for (_,s) in chrom1[chr1]]
    genesStrandsC2 = [s for (_,s) in chrom2[chr2]]

    ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2), genesHomologiesHpSign,(genesNoHomologiesInWindowC1,genesNoHomologiesInWindowC2),genesHomologyGroupsInWindow) = genesComputeHomologyInformations(chrom1, chrom2)

    # Search diagonals
    #FIXME : calculate synteny blocks before, on the whole chromosome, not the ROI specified by the user ranges
    listOfDiags = list(myDiags.extractSbsInPairCompGenomes(chrom1,
                                                           chrom2,
                                                           ancGenes,
                                                           tandemGapMax=arguments['tandemGapMax'],
                                                           gapMax=arguments['gapMax'],
                                                           consistentSwDType=arguments["consistentSwDType"],
                                                           filterType=filterType,
                                                           minChromLength=arguments["minChromLength"],
                                                           distanceMetric=arguments["distanceMetric"],
                                                           pThreshold=arguments["pThreshold"],
                                                           verbose=arguments['verbose']))
    genesDiagIndices = genesComputeDiagIndices(listOfDiags)

    strArray = drawHomologyMatrixWithSBs.drawHomologyMatrix(
        (range1, range2), (genesStrandsC1, genesStrandsC2),
        (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
        (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
        genesHomologiesHpSign, genesHomologyGroupsInWindow,
        genesDiagIndices,
        outputFileName=arguments["out:ImageName"], maxWidth=100, maxHeight=100 )

else:
    #chromosomes are shown as lists of tbs
    g1aID = myMapping.labelWithAncGeneID(genome1, ancGenes)
    g2aID = myMapping.labelWithAncGeneID(genome2, ancGenes)
    ((g1f, gf2gaID1, (nCL1, nGL1)),
     (g2f, gf2gaID2, (nCL2, nGL2))) =\
        myDiags.filter2D(g1aID, g2aID,
                         filterType,
                         arguments["minChromLength"],
                         keepOriginal=True)
    (g1tb, gtb2gf1, nGTD1) =\
        myMapping.remapRewriteInTb(g1f,
                                   tandemGapMax=arguments['tandemGapMax'],
                                   mOld=None)
    (g2tb, gtb2gf2, nGTD2) =\
        myMapping.remapRewriteInTb(g2f,
                                   tandemGapMax=arguments['tandemGapMax'],
                                   mOld=None)
    gtb2gaID1 = {}
    for c in gtb2gf1:
        # see Mapping class addition
        gtb2gaID1[c] = gtb2gf1[c] + gf2gaID1[c]
    gtb2gaID2 = {}
    for c in gtb2gf2:
        # see Mapping class addition
        gtb2gaID2[c] = gtb2gf2[c] + gf2gaID2[c]

    ##print >> sys.stderr, "List of (chromosomes, length in tbs) of Genome 1, for chr of size > %s " % thresholdChr
    ##for (chr1,len1) in [(key1, len(chr1)) for (key1,chr1) in g1tb.items() if len(chr1) > thresholdChr]:
    ##       print >> sys.stderr, "chr %s has %s tbs" % (chr1,len1)
    ##print >> sys.stderr, "List of (chromosomes, length in tbs) of Genome 2, for chr of size > %s " % thresholdChr
    ##for (chr2,len2) in [(key2, len(chr2)) for (key2,chr2) in g2tb.items() if len(chr2) > thresholdChr]:
    ##       print >> sys.stderr, "chr %s has %s tbs" % (chr2,len2)

    ## Inverse dictionnaries for conversion from the gene base to the tb base
    #g12g1filt = dict([(chr1,dict([(i1,i1filt) for (i1filt, i1) in gf2gaID1[chr1].items()])) for chr1 in gf2gaID1])
    #g22g2filt = dict([(chr2,dict([(i2,i2filt) for (i2filt, i2) in gf2gaID2[chr2].items()])) for chr2 in gf2gaID2])
    #g1filt2g1tb = collections.defaultdict(lambda:collections.defaultdict(int))
    #g2filt2g2tb = collections.defaultdict(lambda:collections.defaultdict(int))
    #for chr1 in gtb2gf1:
    #    for (i1tb, i1sfilt) in enumerate(gtb2gf1[chr1]):
    #        for i1filt in i1sfilt:
    #            g1filt2g1tb[chr1][i1filt] = i1tb
    #g2filt2g2tb = collections.defaultdict(lambda:collections.defaultdict(int))
    #for chr2 in gtb2gf2:
    #    for (i2tb, i2sfilt) in enumerate(gtb2gf2[chr2]):
    #        for i2filt in i2sfilt:
    #            g2filt2g2tb[chr2][i2filt] = i2tb

    if not arguments['convertGenicToTbCoordinates']:
        (chr1,range1) = parseChrRange(arguments["chr1:deb1-fin1"], g1tb)
        (chr2,range2) = parseChrRange(arguments["chr2:deb2-fin2"], g2tb)
    else:

        #if the user uses genic coordinates and wants to view results in the tb coordinates
        #inverse again dictionnaries
        #g12g1tb = collections.defaultdict(lambda:collections.defaultdict(int))
        #for chr1 in g1filt2g1tb:
        #    i1_old=None
        #    i1_old_tmp=0
        #    for i1,_ in enumerate(g1aID[chr1]):
        #        if i1 in g12g1filt[chr1]:
        #            if i1_old == None:
        #                for i11 in range(i1_old_tmp):
        #                    g12g1tb[chr1][i11]=g1filt2g1tb[chr1][g12g1filt[chr1][i1]]
        #            g12g1tb[chr1][i1]=g1filt2g1tb[chr1][g12g1filt[chr1][i1]]
        #            i1_old = i1
        #        else:
        #            if i1_old == None:
        #                i1_old_tmp+=1
        #            else:
        #                g12g1tb[chr1][i1]=g1filt2g1tb[chr1][g12g1filt[chr1][i1_old]]
        #g22g2tb = collections.defaultdict(lambda:collections.defaultdict(int))
        #for chr2 in g2filt2g2tb:
        #    i2_old=None
        #    i2_old_tmp=0
        #    for i2,_ in enumerate(g2aID[chr2]):
        #        if i2 in g22g2filt[chr2]:
        #            if i2_old == None:
        #                for i22 in range(i2_old_tmp):
        #                    g22g2tb[chr2][i22]=g2filt2g2tb[chr2][g22g2filt[chr2][i2]]
        #            g22g2tb[chr2][i2]=g2filt2g2tb[chr2][g22g2filt[chr2][i2]]
        #            i2_old = i2
        #        else:
        #            if i2_old == None:
        #                i2_old_tmp+=1
        #            else:
        #                g22g2tb[chr2][i2]=g2filt2g2tb[chr2][g22g2filt[chr2][i2_old]]
        #(chr1,range1) = parseChrRange(arguments["chr1:deb1-fin1"], g1tb, g2gtb=g12g1tb)
        #(chr2,range2) = parseChrRange(arguments["chr2:deb2-fin2"], g2tb, g2gtb=g22g2tb)
        (chr1,range1) = parseChrRange(arguments["chr1:deb1-fin1"], g1tb, g2gtb=gtb2gaID1.old)
        (chr2,range2) = parseChrRange(arguments["chr2:deb2-fin2"], g2tb, g2gtb=gtb2gaID2.old)


    chrom1_tb ={}
    chrom2_tb ={}
    chrom1_tb[chr1] = g1tb[chr1][range1[0]:range1[1]]
    chrom2_tb[chr2] = g2tb[chr2][range2[0]:range2[1]]

    #Focus on the chromosome of the window, just give simple name to the chromosome of interest
    chrom1 = genome1[chr1]
    c1_aID = g1aID[chr1]
    c1_aID_filt = g1f[chr1]
    Cf2CaID1 = gf2gaID1[chr1]
    Ctb2Cf1 = gtb2gf1[chr1]
    chrom2 = genome2[chr2]
    c2_aID = g2aID[chr2]
    c2_aID_filt = g2f[chr2]
    Cf2CaID2 = gf2gaID2[chr2]
    Ctb2Cf2 = gtb2gf2[chr2]
    #chrom12chrom1tb = g12g1tb[chr1] #useless
    #chrom22chrom2tb = g22g2tb[chr2] #useless

    ###
    # Build TbNumberOfGenesInEachTbC1 : [ 4,5,1,1,6,2, ...] number og genes in each TB of C1
    ###
    TbNumberOfGenesInEachTbC1 = [len(Ctb2Cf1[i1_tb]) for i1_tb in range(range1[0],range1[1])]
    TbNumberOfGenesInEachTbC2 = [len(Ctb2Cf2[i2_tb]) for i2_tb in range(range2[0],range2[1])]

    ###
    # Build TBStrands
    ###
    TbStrandsC1 = [s for (_,s) in chrom1_tb[chr1]]
    TbStrandsC2 = [s for (_,s) in chrom2_tb[chr2]]

    ###
    # Build rangeXTB
    ###
    (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow) =\
        TbComputeHomologyInformations(chrom1_tb[chr1], chrom2_tb[chr2])

    # compute the recommended gapMax parameter
    #########################################
    N12s, N12_g = myDiags.numberOfHomologies(chrom1_tb, chrom2_tb)
    (p_hpSign,p_hpSign_g,(sTBG1, sTBG1_g),(sTBG2, sTBG2_g)) = myProbas.statsHpSign(chrom1_tb,chrom1_tb)
    p_hpSign={(chr1,chr2):p_hpSign[(1,1)]}#FIXME
    print >> sys.stderr, "genome 1 tb orientation proba = {+1:%s,-1:%s,None:%s} (stats are also calculated for each chromosome)" % (sTBG1_g[+1], sTBG1_g[-1], sTBG1_g[None])
    print >> sys.stderr, "genome 2 tb orientation proba = {+1=%s,-1:%s,None:%s} (stats are also calculated for each chromosome)" % (sTBG2_g[+1], sTBG2_g[-1], sTBG2_g[None])
    print >> sys.stderr, "hp sign proba in the 'global' mhp = {+1:%s,-1:%s,None:%s) (probabilities are for pairwise mhp)" % (p_hpSign_g[+1], p_hpSign_g[-1], p_hpSign_g[None])
    N1_g=sum([len(chrom1_tb[c1]) for c1 in chrom1_tb])
    N2_g=sum([len(chrom2_tb[c2]) for c2 in chrom2_tb])
    m=2
    #FIXME : calculate adviced gapMax before, on the whole chromosome, not the ROI specified by the user ranges
    gap = myDiags.recommendedGap(arguments["nbHpsRecommendedGap"],
                                 arguments["targetProbaRecommendedGap"],
                                 N12_g,
                                 N1_g,
                                 N2_g,
                                 p_hpSign=p_hpSign_g,
                                 verbose=arguments['verbose'])
    print >> sys.stderr, "recommended gapMax = %s tbs" % gap
    if arguments['gapMax'] == None:
        gapMax = gap
    else:
        gapMax = arguments['gapMax']
    print >> sys.stderr, "used gapMax = %s" % gapMax

    # Search diagonals
    listOfDiags = myDiags.extractSbsInPairCompChr(chr1,
                                                  chr2,
                                                  chrom1_tb[chr1],
                                                  chrom2_tb[chr2],
                                                  consistentSwDType=arguments['consistentSwDType'],
                                                  gapMax=arguments['gapMax'],
                                                  distanceMetric=arguments["distanceMetric"],
                                                  verbose=arguments['verbose'])
    # statistical validation (SV)
    #FIXME : calculate p-values before, on the whole chromosome, not the ROI specified by the user ranges
    print >> sys.stderr, "Warning the p-value calculation is performed on the ROI defined by the user ranges"
    listOfDiags = list(myDiags.statisticalValidation(listOfDiags,
                                                           chrom1_tb,
                                                           chrom2_tb,
                                                           N12s,
                                                           p_hpSign,
                                                           pThreshold = arguments['pThreshold'],
                                                           NbOfHomologiesThreshold=50,
                                                           verbose=arguments['verbose']))
    TbDiagIndices = TbComputeDiagIndices(listOfDiags)

    strArray =\
        drawHomologyMatrixWithSBs.drawHomologyMatrix((range1, range2),
                                                     (TbStrandsC1, TbStrandsC2),
                                                     ([],[]),
                                                     (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2),
                                                     TbHpSign,
                                                     TbHomologyGroupsInWindow,
                                                     TbDiagIndices,
                                                     outputFileName=arguments["out:ImageName"],
                                                     maxWidth=100, maxHeight=100,
                                                     symbolsInGenes=(TbNumberOfGenesInEachTbC1, TbNumberOfGenesInEachTbC2))

#copy the css style sheet
dirNameImage = os.path.dirname(arguments["out:ImageName"])
dirNameImage = dirNameImage if dirNameImage != "" else "."
print >> sys.stderr, "cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage)
os.system("cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage))

# write a simple file with all diagonals into output file
f=open(arguments['out:SyntenyBlocks'],'w')
print >> f, "Mode : %s" %  'Genic scale' if arguments['mode:chromosomesRewrittenInTbs'] == False else 'Tandem Blocks scale'
print >> f, "chromosome %s de %s\t%s\t%s\tchromosome %s de %s\t%s\t%s\t%s" % (chr1, genome1Name, 'beginC1', 'endC1', chr2, genome2Name, 'beginC2', 'endC2','length in ancestral genes')
print >> f, "c1\tbeg1\tend1\tc2\tbeg2\tend2\thps\tpVal"
listOfDiags = [l for l in listOfDiags]
listOfDiags = sorted(list(listOfDiags), key=lambda x:len(x[2]), reverse=True)
for diag in listOfDiags:
    ((c1,l1), (c2,l2), la, pVal) = diag
    min_l2 = min(l2[0],l2[-1])
    max_l2 = max(l2[0],l2[-1])
    print >> f, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (c1,l1[0],l1[-1],c2,min_l2,max_l2,len(la),pVal)
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
if not arguments['mode:chromosomesRewrittenInTbs']:
    title =\
        "%s, tandemGapMax=%s tbs, gapMax=%s%s, %s sbs" %\
        ('MH',
         arguments['tandemGapMax'],
         arguments['gapMax'],
         arguments['distanceMetric'],
         len(listOfDiags))
else:
    title =\
        "%s, tandemGapMax=%s tbs, gapMax=%s%s, %s sbs" %\
        ('MHP',
         arguments['tandemGapMax'],
         arguments['gapMax'],
         arguments['distanceMetric'],
         len(listOfDiags))\

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

#Find diags with the more paralogs
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
    l1=diag[0][1]
    l2=diag[1][1]
    la=diag[2]
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
def chooseChrsAndRanges(genome1, genome2, ancGenes, distanceMetric = 'DPD'):
    while True:
        try:
            (chr1,range1) = parseChrRange(raw_input("chr1:deb1-fin1 = "), genome1)
            break
        except ValueError:
            print >> sys.stderr, "You need to write somtehing as chr1:deb1-fin1 with chr1 a chr of G1 and deb1 and fin1 indices of the first and last gene (indices start at 1)"

    while True:
        try:
            (chr2,range2) = parseChrRange(raw_input("chr2:deb2-fin2 = "), genome2)
            break
        except ValueError:
            print >> sys.stderr, "You need to write somtehing as chr2:deb2-fin2 with chr2 a chr of G2 and deb2 and fin2 indices of the first and last gene (indices start at 1)"

    chrom1 ={}
    chrom2 ={}
    chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
    chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]
    listOfDiags = myDiags.extractSbsInPairCompGenomes(chrom1, chrom2, ancGenes, gapMax=arguments["gapMax"], consistentSwDType=arguments["consistentSwDType"], filterType=filterType, minChromLength=arguments["minChromLength"], distanceMetric=arguments["distanceMetric"], verbose=arguments['verbose'])
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
        return chooseChrsAndRanges(genome1, genome2, ancGenes)
    else:
        return (chrom1, chrom2, range1, range2)
