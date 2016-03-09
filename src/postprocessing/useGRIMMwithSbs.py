#!/usr/bin/python
# -*- coding: utf-8 -*-

# This script wraps the GRIMM v2.01 software
# http://grimm.ucsd.edu/GRIMM/

# # Install GRIMM
# PATH_GRIMM_PARENT=/home/${USER}/Libs
# mkdir -p ${PATH_GRIMM}
# cd ${PATH_GRIMM}
# wget http://grimm.ucsd.edu/DIST/GRIMM-2.01.zip
# unzip GRIMM-2.01.zip
# PATH_GRIMM=${PATH_GRIMM_PARENT}/GRIMM-2.01
# cd ${PATH_GRIMM}
# # Compile
# make
# # If compilation does not work, comment these lines in ./MakeFile:
# ---------------------------------------
# ifeq ($(strip $(CC)), gcc)
# ifeq ($(OS), Linux)
# CFLAGS := -march=pentiumpro $(CFLAGS)
# endif
# ifeq ($(OS), SunOS)
# CFLAGS := -mv8 $(CFLAGS)
# endif
# endif
# ----------------------------------------
# echo "export PATH=\"\$PATH:${PATH_GRIMM}\"" >> ~/.bashrc
# cd
#

# The aim is to use PhylDiag's conserved segments between human and mouse in the X chromosome
# and use these conserved segments to retrieve reversal scenarios inferred from GRIMM

import collections
import itertools
from numpy import random
import sys
import os
import subprocess

import math

from utils import myTools, myDiags, myLightGenomes, myGenomes, myGenomesDrawer, mySvgDrawer, myIntervals
import bidict
from utils. mySvgDrawer import Point as Point
import libs.myEvents as myEvents
import libs.myBreakpointsAnalyser

# TODO : compute the reuse statistic r (Sankoff 2004)
# 1 < r = 2d/b < 2 avec b=b'-c
# b' le nombre de block sur les c chromosomes
# c le nb de chromosomes used as input (b' sbs on it)
#
GENOMENAME = 0

def projectSbsOnGenome(listOfSbsInCompC1C2, cx):
    assert cx in [1, 2]
    # give orientations to each sb
    if cx == 1:
        # orientations of blocks are all set by default to +1 in the reference chromosome, except when the sb has a None diagType
        listSbsInCx = [(id, sb, (+1 if sb.dt is not None else None)) for (sb, id) in listOfSbsInCompC1C2]
    else:
        assert cx == 2
        listSbsInCx = []
        for (i, (sb, id)) in enumerate(listOfSbsInCompC1C2):
            if sb.dt == '/':
                s = +1
            elif sb.dt == '\\':
                s = -1
            else:
                assert sb.dt == None
                s = None
            listSbsInCx.append((id, sb, s))
    return listSbsInCx

def rewriteGenomesIntoSbs(sbsInPairComp, checkNoOverlap=True, inputWithIds=True, reverseRefGenome=False):

    if inputWithIds == False:
        # give ids to each sbs
        assert isinstance(sbsInPairComp, myTools.Dict2d)
        sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
        for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
            assert isinstance(sb, myDiags.SyntenyBlock)
            sbsInPairCompWithIds.addToLocation((c1, c2), sb)
    else:
        assert isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists)
        sbsInPairCompWithIds = sbsInPairComp

    genome1 = myLightGenomes.LightGenome()
    genome2 = myLightGenomes.LightGenome()
    # [for (id, (c1, c2), s) in sbsInPairCompWithIds.iterByOrderedIds()]
    listOfSortedSbsInC1 = collections.defaultdict(list)
    listOfSortedSbsInC2 = collections.defaultdict(list)
    for (c1, c2) in [(c1, c2) for c1 in sbsInPairCompWithIds for c2 in sbsInPairCompWithIds[c1]]:
        listOfSortedSbsInC1[c1] += projectSbsOnGenome(sbsInPairCompWithIds.getItemsAndIdsByLocation((c1, c2)), 1 if not reverseRefGenome else 2)
        # FIXME: sbs of ...InC1 are not good
        listOfSortedSbsInC2[c2] += projectSbsOnGenome(sbsInPairCompWithIds.getItemsAndIdsByLocation((c1, c2)), 2 if not reverseRefGenome else 1)
        assert len(listOfSortedSbsInC2[c2]) > 0

    for (c1, c2) in [(c1, c2) for c1 in sbsInPairCompWithIds for c2 in sbsInPairCompWithIds[c1]]:
        # sort sbs by starting position on each genome
        # listOfSortedSbsInC1[c1] = [..., (id, sb, s), ...]
        listOfSortedSbsInC1[c1].sort(key=lambda x: x[1].minOnG(1))
        listOfSortedSbsInC2[c2].sort(key=lambda x: x[1].minOnG(2))

        if checkNoOverlap:
            for ((id1, sb1, s1), (id2, sb2, s2)) in myTools.myIterator.slidingTuple(listOfSortedSbsInC1[c1]):
                assert sb1.minOnG(1) < sb2.minOnG(1)
                if sb2.minOnG(1) < sb1.maxOnG(1):
                    print >> sys.stderr, 'Warning, overlap of sbs: ' + "sb2.minOnG(1)=%s < sb1.maxOnG(1)=%s" % (sb2.minOnG(1), sb1.maxOnG(1))
                    print >> sys.stderr, 'sb1.l1=%s' % sb1.l1
                    print >> sys.stderr, 'sb2.l1=%s' % sb2.l1
                # assert sb1.maxOnG(1) < sb2.minOnG(1), "%s < %s" % (sb1.maxOnG(1), sb2.minOnG(1))
            for ((id1, sb1, s1), (id2, sb2, s2)) in myTools.myIterator.slidingTuple(listOfSortedSbsInC2[c2]):
                assert sb1.minOnG(2) < sb2.minOnG(2)
                #assert sb1.maxOnG(2) < sb2.minOnG(2), "%s < %s" % (sb1.maxOnG(2), sb2.minOnG(2))
                if sb2.minOnG(2) < sb1.maxOnG(2):
                    print >> sys.stderr, 'Warning, overlap of sbs: ' + "sb2.minOnG(2)=%s < sb1.maxOnG(2)=%s" % (sb2.minOnG(2), sb1.maxOnG(2))
                    print >> sys.stderr, 'sb1.l1=%s' % sb1.l2
                    print >> sys.stderr, 'sb2.l1=%s' % sb2.l2

        genome1[c1] = [myLightGenomes.OGene(str(id), s) for (id, sb, s) in listOfSortedSbsInC1[c1]]
        genome2[c2] = [myLightGenomes.OGene(str(id), s) for (id, sb, s) in listOfSortedSbsInC2[c2]]
    return (genome1, genome2)

def printGenomeInto_GRIMM_format(genome, stream=sys.stdout):
    assert isinstance(stream, file)
    global GENOMENAME
    if genome.name is None:
        GENOMENAME += 1
        genome.name = str(GENOMENAME)
    print >> stream, ">%s" % genome.name
    for (c, chrom) in genome.iteritems():
        print >> stream, "# Chromosome %s" % c
        chrStr = ''
        for g in chrom:
            if g.s == +1:
                s = ''
            elif g.s == -1:
                s = '-'
            else:
                assert g.s is None
                s = None
            if s is not None:
                chrStr += "%s%s " % (s, g.n)
            else:
                chrStr += "[ %s ] " % g.n
        print >> stream, chrStr + '$'

def main():
    __doc__ = \
    """use synteny blocks of PhylDiag to find the minimum scenario of rearrangements from GRIMM
    """

    arguments = myTools.checkArgs([('syntenyBlocks', file), ('genome1', file), ('genome2', file)],
                                  [('outPath', str, './'),
                                   ('pathToGRIMM', str, '/home/jlucas/Libs/GRIMM-2.01/grimm'),
                                   ('renameSbsByOrder', bool, True),
                                   ('removeSbsWithAnUnknownOrientation', bool, True),
                                   ('selectedChrom', str, 'X'),
                                   ('unitOfGenomes', str, 'Mb'),
                                   ('removeMonoGenicConservedSegments', bool, False),
                                   # show the scenario from genome2 to genome 1 if True
                                   ('invertSpeciesOrder', bool, False)],
                                  __doc__)
    # TODO: unit of genomes = gene for simu (first and last genomes)
    # FIXME: lengths of genes for simu

    assert arguments['unitOfGenomes'] in {'gene', 'Mb'}
    unitOfGenomes = arguments['unitOfGenomes']

    genome1 = myLightGenomes.LightGenome(arguments['genome1'], withDict=True)
    genome2 = myLightGenomes.LightGenome(arguments['genome2'], withDict=True)
    genome1H = myGenomes.Genome(arguments['genome1'])
    genome2H = myGenomes.Genome(arguments['genome2'])
    sbsInPairCompWithIds = myDiags.parseSbsFile(arguments['syntenyBlocks'],
                                                genome1=genome1,
                                                genome2=genome2,
                                                withIds=True)
    if arguments['invertSpeciesOrder']:
        (genome1, genome2) = (genome2, genome1)
        (genome1H, genome2H) = (genome2H, genome1H)
        for (id, (k1, k2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
            (sb.l1, sb.l2) = (sb.l2, sb.l1)

    sChr = arguments['selectedChrom']
    # Only consider sbs in the selected chromosome
    new_sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    for (id, (k1, k2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        if k1 == sChr and k2 == sChr:
            assert isinstance(sb, myDiags.SyntenyBlock)
            if arguments['removeMonoGenicConservedSegments'] and len(sb.la) == 1:
                continue
            new_sbsInPairCompWithIds.addToLocationWithId((k1, k2), sb, id)
    sbsInPairCompWithIds = new_sbsInPairCompWithIds

    # rename : id=1 (longest sb), .....
    if arguments['renameSbsByOrder']:
        new_sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
        newId = 1
        for (id, (k1, k2), sb) in sorted([(id, (k1, k2), sb) for (id, (k1, k2), sb) in sbsInPairCompWithIds.iterByOrderedIds()], key=lambda x: len(x[2].la), reverse=True):
                new_sbsInPairCompWithIds.addToLocationWithId((k1, k2), sb, newId)
                newId += 1
        sbsInPairCompWithIds = new_sbsInPairCompWithIds

    # record these synteny blocks
    with open(arguments['outPath'] + '/syntenyBlocksOnSelectedChr.txt', 'w') as f:
        myDiags.printSbsFile(sbsInPairCompWithIds, genome1, genome2, stream=f)

    # get the lengths of sbs
    sbId2MeanLengthU = myDiags.getSbsMeanLengths(genome1H, genome2H, sbsInPairCompWithIds, lengthUnit=unitOfGenomes)

    # FIXME does not work for magSimus genomes beg=i, end=i+1
    sbIdCoordsInUInGenome1 = collections.defaultdict(lambda: collections.defaultdict(tuple))
    sbIdCoordsInUInGenome2 = collections.defaultdict(lambda: collections.defaultdict(tuple))
    def findSbCoordsInUOnGenome(sb, c, genomeH, rankGenome, unit=unitOfGenomes, defaultHalfIntergeneLength=None, defaultGeneLength=None):
        assert rankGenome in {1, 2}
        assert isinstance(sb, myDiags.SyntenyBlock)
        assert unitOfGenomes in {'Mb', 'gene'}
        l = sb.l1 if rankGenome == 1 else sb.l2
        if unitOfGenomes == 'gene':
            assert isinstance(defaultHalfIntergeneLength, int) or isinstance(defaultHalfIntergeneLength, float)
            assert isinstance(defaultGeneLength, int) or isinstance(defaultHalfIntergeneLength, float)
            minOnG = min([idxG for idxGs in l for idxG in idxGs])
            maxOnG = max([idxG for idxGs in l for idxG in idxGs])
            intergeneLengthsBefore = (2 * minOnG + 1) * defaultHalfIntergeneLength
            interGeneLengthInside = ((maxOnG - minOnG) * 2 + 1 ) * defaultHalfIntergeneLength
            geneLengthsBefore = minOnG * defaultGeneLength
            geneLengthsInside = (maxOnG - minOnG) * defaultGeneLength
            lengthBefore = intergeneLengthsBefore + geneLengthsBefore
            lengthInside = interGeneLengthInside + geneLengthsInside
            res = (lengthBefore, lengthBefore + lengthInside)
        else:
            assert unitOfGenomes == 'Mb'
            genesInL = [genomeH.lstGenes[c][idxG] for idxGs in l for idxG in idxGs]
            allCoordsOfGenesInL = []
            for gene in genesInL:
                assert isinstance(gene, myGenomes.Gene)
                assert (gene.strand in {+1, None} and gene.beginning < gene.end) or (gene.strand == -1 and gene.end< gene.beginning)
                assert isinstance(gene, myGenomes.Gene)
                allCoordsOfGenesInL.append(gene.beginning)
                allCoordsOfGenesInL.append(gene.end)
            (minOnG, maxOnG) = (min(allCoordsOfGenesInL), max(allCoordsOfGenesInL))
            res = (minOnG / 1000000.0, maxOnG / 1000000.0)
        return res


    halfIntergeneLength = 1.5
    defaultGeneLength = 1.0
    defaultHalfIntergene = 0.5
    for (id, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        sbIdCoordsInUInGenome1[c1][id] = findSbCoordsInUOnGenome(sb, myGenomes.commonChrName(c1), genome1H, 1, unit=unitOfGenomes,
                                                                 defaultHalfIntergeneLength=defaultHalfIntergene, defaultGeneLength=defaultGeneLength)
        sbIdCoordsInUInGenome2[c2][id] = findSbCoordsInUOnGenome(sb, myGenomes.commonChrName(c2), genome2H, 2, unit=unitOfGenomes,
                                                                 defaultHalfIntergeneLength=defaultHalfIntergene, defaultGeneLength=defaultGeneLength)

    (genome1InSbs, genome2InSbs) = rewriteGenomesIntoSbs(sbsInPairCompWithIds, reverseRefGenome=arguments['invertSpeciesOrder'])
    genome1InSbs.name = genome1.name
    genome2InSbs.name = genome2.name
    genome1InSbs.sort(byName=True)
    genome2InSbs.sort(byName=True)
    # # print these two representation of the genomes
    # pathToSyntenyBlocksForGrimm = arguments['outPath'] + '/syntenyBlocksForGRIMM_oldNames.txt'
    # with open(pathToSyntenyBlocksForGrimm, 'w') as f:
    #     printGenomeInto_GRIMM_format(genome1InSbs, stream=f)
    #     printGenomeInto_GRIMM_format(genome2InSbs, stream=f)

    # reduce genomes to only chromosomes X
    genome1OnlySchr = myLightGenomes.LightGenome()
    genome2OnlySchr = myLightGenomes.LightGenome()
    genome1OnlySchr.name = genome1.name
    genome2OnlySchr.name = genome2.name
    sChr = arguments['selectedChrom']
    genome1OnlySchr[sChr] = genome1InSbs[sChr]
    genome2OnlySchr[sChr] = genome2InSbs[sChr]
    assert len(genome1OnlySchr[sChr]) > 0
    # FIXME is it usefull ?
    geneNamesInGenomeAndNotTheOther = genome1OnlySchr.getGeneNames() ^ genome2OnlySchr.getGeneNames()
    genome1OnlySchr.removeGenes(geneNamesInGenomeAndNotTheOther)
    genome2OnlySchr.removeGenes(geneNamesInGenomeAndNotTheOther)
    assert len(genome1OnlySchr[sChr]) == len(genome2OnlySchr[sChr])

    # renaming should be done because of internal limitations of GRIMM (that do not want gene names to be superior than
    # the total number of genes.
    geneNamesInBothGenomes = genome1OnlySchr.getGeneNames()
    assert geneNamesInBothGenomes == genome2OnlySchr.getGeneNames()
    oldGn2newGn = bidict.bidict()
    if max([int(gn) for (gn, gs) in genome1OnlySchr[sChr]]) > len(genome1OnlySchr[sChr]):
        newGeneName = 1
        for gn in geneNamesInBothGenomes:
            oldGn2newGn[gn] = str(newGeneName)
            newGeneName += 1
        for i, g in enumerate(genome1OnlySchr[sChr]):
            genome1OnlySchr[sChr][i] = myLightGenomes.OGene(oldGn2newGn[g.n], g.s)
        for i, g in enumerate(genome2OnlySchr[sChr]):
            genome2OnlySchr[sChr][i] = myLightGenomes.OGene(oldGn2newGn[g.n], g.s)
    else:
        for gn in geneNamesInBothGenomes:
            oldGn2newGn[gn] = gn

    # check that there are no conserved adjacencies of cs
    setOfOA1 = myIntervals.analyseGenomeIntoAdjacencies(genome1OnlySchr)
    setOfOA2 = myIntervals.analyseGenomeIntoAdjacencies(genome2OnlySchr)
    setOfConservedOAdjs = setOfOA1 & setOfOA2
    if len(setOfConservedOAdjs) > 0:
        print >> sys.stderr, 'WARNING !!!!!!!!!!'
        print >> sys.stderr, 'The is at least one conserved adjaceny of cs between the two genome before computing the parcimonious scenario with GRIMM'
        print >> sys.stderr, 'the conserved adjs are : %s' % [(oldGn2newGn.inv[oadj.n], oadj.s) for oadj in setOfConservedOAdjs]
        print >> sys.stderr, 'WARNING !!!!!!!!!!'
    pathToSyntenyBlocksForGrimm = arguments['outPath'] + '/syntenyBlocksForGRIMM_specificNames.txt'
    with open(pathToSyntenyBlocksForGrimm, 'w') as f:
        printGenomeInto_GRIMM_format(genome1OnlySchr, stream=f)
        printGenomeInto_GRIMM_format(genome2OnlySchr, stream=f)
    # lengths in bp of all conserved segments
    lengthOfCS = {}

    # wrap the grimm executable
    pathToGrimm = arguments['pathToGRIMM']
    if os.path.isfile(pathToGrimm):
        options = ' -v -L'
        bashComand = pathToGrimm + options +' -f ' + pathToSyntenyBlocksForGrimm
        subProcStdin = None  #  subprocess.PIPE, if you want it
        subProcStderr = None  #  subprocess.PIPE, if you want it
        proc = subprocess.Popen(bashComand, shell=True, stdin=subProcStdin, stdout=subprocess.PIPE, stderr=subProcStderr)
        # returns stdout and stderr
        grimmStdOut, subProcStderr = proc.communicate(subProcStdin)
        # proc.returncode, might be interesting too
        # stdout is a str
        # Use this line to output the output of GRIMM
        with open(arguments['outPath'] + '/grimmOutputFile.txt', 'w') as grimmOutputFile:
            print >> grimmOutputFile, grimmStdOut

        approx_equal = lambda a, b, t: abs(a - b) < t
        stepWasOnPrevLine = False
        # print >> sys.stderr, [oldGn2newGn.inv[gn] for (gn, s) in genome1OnlySchr[sChr]]
        # print >> sys.stderr, sbId2Length.keys()
        totalLength = sum([sbId2MeanLengthU[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr]])
        for (n, s) in genome1OnlySchr[sChr]:
            assert s in{+1, -1, None}

        seevolutionInputFile = open(arguments['outPath'] + '/seevolutionInputFile.txt', 'w')
        # input xml file for seevolution
        print >> seevolutionInputFile, '<?xml version="1.0" encoding="ISO-8859-1"?>'
        print >> seevolutionInputFile, '<genomemutations id="1">'
        print >> seevolutionInputFile, '<organism id="Organism">'
        # the small +1 is here to avoid bound error in seevolution
        # only the circular version works
        print >> seevolutionInputFile, '<chromosome id="%s" circular="true" length="%s"/>' % (sChr, int(round(totalLength, 0)))
        print >> seevolutionInputFile, '<speciation name="ancestor">'

        genome1OnlySchrWithOldNames = myLightGenomes.LightGenome(genome1OnlySchr)
        genome1OnlySchrWithOldNames[sChr] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome1OnlySchr[sChr]]
        genome2OnlySchrWithOldNames = myLightGenomes.LightGenome(genome2OnlySchr)
        genome2OnlySchrWithOldNames[sChr] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome2OnlySchr[sChr]]

        sbCoordsInUChr1 = [sbIdCoordsInUInGenome1[sChr][int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr]]
        sbCoordsInUChr2 = [sbIdCoordsInUInGenome2[sChr][int(gn)] for (gn, s) in genome2OnlySchrWithOldNames[sChr]]

        widthGene = 100.0 / 100.0
        desiredLength = 100.0

        colorGenerator = myGenomesDrawer.LevelIdxGenerator()
        name2color = {}
        for gn in [g.n for g in genome1OnlySchrWithOldNames[sChr]]:
            assert gn not in name2color
            name2color[gn] = colorGenerator.getLevel()

        symbolsInGenes = [g.n for g in genome1OnlySchrWithOldNames[sChr]]
        (items_genomeWithSbCoords1, lengthFactor) = myGenomesDrawer.svgItemsChromosome(genome1OnlySchrWithOldNames[sChr],
                                                                                       colorsGenerator=None,
                                                                                       genesCoordinates=sbCoordsInUChr1,
                                                                                       lengthGenes=None,
                                                                                       halfIntergeneLengths=None,
                                                                                       desiredLength=desiredLength,
                                                                                       symbolsInGenes=symbolsInGenes,
                                                                                       geneWidth=widthGene,
                                                                                       name2color=name2color)
        cptStep = 1
        scenarioItems = []
        radiusBR = 0.65
        BRopacity = 0.5
        interLineFactor = 3.0
        yOrigin = 1 * interLineFactor
        sceneItems = mySvgDrawer.translateItems(items_genomeWithSbCoords1, (0, yOrigin))

        bra = libs.myBreakpointsAnalyser.BreakpointsAnalyser(genome1OnlySchr)
        nbBreakpoints = 0
        nbBreakpointReuse1 = 0
        nbBreakpointReuse2 = 0
        nbBreakpointsAtChromExtremity = 0
        setRemnantBreakpointFlanks = set([])
        setRemnantBreakpointFlanks2 = set([])
        dictBreakpointReusePerFlank = collections.defaultdict(int)
        setGeneExtremitiesWithBreakpointReuse = set()

        for line in grimmStdOut.split('\n'):
            print >> sys.stderr, line
            if 'Step' in line:
                stepWasOnPrevLine = True
                if '(Source)' in line:
                    symbolsInGenes = [g.n for g in genome1OnlySchrWithOldNames[sChr]]
                    lengthGenes = [sbId2MeanLengthU[int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr]]
                    (scenarioItems, factorLength) = myGenomesDrawer.svgItemsChromosome(genome1OnlySchrWithOldNames[sChr],
                                                           lengthGenes=lengthGenes,
                                                           halfIntergeneLengths=halfIntergeneLength,
                                                           desiredLength=desiredLength,
                                                           symbolsInGenes=symbolsInGenes,
                                                           geneWidth=widthGene,
                                                           name2color=name2color)
                    scenarioItems = mySvgDrawer.translateItems(scenarioItems, (0, interLineFactor * cptStep))
                    pass
                else:
                    linesBracketsL = line.split('[')
                    linesGBL = line.split('gene ')
                    # leftGn: leftmost gene on the left of the inverted segment
                    (leftGn, rightGn) = (linesBracketsL[1].split(']')[0],  linesBracketsL[2].split(']')[0])
                    # rightGn: rightmost gene on the right of the inverted segment
                    (leftGn, rightGn) = (leftGn.replace('-', ''), rightGn.replace('-', ''))
                    # to get the former id:
                    # leftGn_old = oldGn2newGn.inv[leftGn]
                    # indexes start at 0
                    (leftGidx, rightGidx) = (linesGBL[1].split(' [')[0],  linesGBL[2].split(' [')[0])
                    (leftGidx, rightGidx) = (int(leftGidx), int(rightGidx))
                    assert leftGidx <= rightGidx
                    # print >> sys.stderr, '(leftGn, rightGn)= (%s, %s)' % (leftGn, rightGn)
                    assert genome1OnlySchr[sChr][leftGidx].n == leftGn, '%s == %s' % (genome1OnlySchr[sChr][leftGidx].n, leftGn)
                    assert genome1OnlySchr[sChr][rightGidx].n == rightGn, '%s == %s' % (genome1OnlySchr[sChr][rightGidx].n, rightGn)

                    leftLengthU = sum([sbId2MeanLengthU[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr][:leftGidx]])
                    invLengthU = sum([sbId2MeanLengthU[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr][leftGidx:rightGidx+1]])
                    rightLengthU = sum([sbId2MeanLengthU[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr][rightGidx+1:]])
                    assert approx_equal(leftLengthU + invLengthU + rightLengthU, totalLength, 0.01)

                    # because of seevolution
                    if leftGidx == rightGidx:
                        invLengthU = 0
                    print >> seevolutionInputFile, "<inversion left=\"%s\" right=\"%s\" chromosome=\"%s\"/>" % (int(round(leftLengthU, 0)),
                                                                                                      int(round(leftLengthU, 0)) + int(round(invLengthU,0)),
                                                                                                      sChr)
                    cptStep += 1

                    # analyse breakpoints before performing them
                    nbBreakpoints += 2
                    assert 0 <= leftGidx < len(genome1OnlySchr[sChr])
                    leftBreakpointFlanks = myIntervals.geneExtremitiesFromIntergeneInChrom(leftGidx, genome1OnlySchr[sChr])
                    reusedLeftFlanksOfBreakpoint = leftBreakpointFlanks & setRemnantBreakpointFlanks
                    if len(reusedLeftFlanksOfBreakpoint) > 0:
                        setGeneExtremitiesWithBreakpointReuse |= reusedLeftFlanksOfBreakpoint
                        nbBreakpointReuse1 += 1
                    assert 0 <= rightGidx < len(genome1OnlySchr[sChr])
                    rightBreakpointFlanks = myIntervals.geneExtremitiesFromIntergeneInChrom(rightGidx + 1, genome1OnlySchr[sChr])
                    reusedRightFlanksOfBreakpoint = rightBreakpointFlanks & setRemnantBreakpointFlanks
                    if len(reusedRightFlanksOfBreakpoint) > 0:
                        setGeneExtremitiesWithBreakpointReuse |= reusedRightFlanksOfBreakpoint
                        nbBreakpointReuse1 += 1
                    setRemnantBreakpointFlanks |= leftBreakpointFlanks | rightBreakpointFlanks

                    # draw the inversion
                    # 0.1 is the halfintergene length
                    leftLengthGraph = factorLength * sum([(sbId2MeanLengthU[int(gn)] + 2 * halfIntergeneLength) for (gn, s) in genome1OnlySchr[sChr][:leftGidx]])
                    invLengthGraph = factorLength * sum([(sbId2MeanLengthU[int(gn)] + 2 * halfIntergeneLength) for (gn, s) in genome1OnlySchr[sChr][leftGidx:rightGidx+1]])
                    rightLengthGraph = leftLengthGraph + invLengthGraph
                    # leftLengthGraph -= halfIntergeneLength / 2.0
                    # rightLengthGraph -= halfIntergeneLength / 2.0
                    # add lines to visualise the inversion
                    scenarioItems.append(mySvgDrawer.Line(Point(leftLengthGraph, interLineFactor * (cptStep-1)), Point(rightLengthGraph, interLineFactor * cptStep), width=0.2))
                    scenarioItems.append(mySvgDrawer.Line(Point(rightLengthGraph, interLineFactor * (cptStep-1)), Point(leftLengthGraph, interLineFactor * cptStep), width=0.2))
                    invLengthInGenes = (rightGidx + 1) - leftGidx
                    # write the length (in genes and in bps) of the inversion on the left
                    if unitOfGenomes == 'Mb':
                        invLengthInMb = sum([sbId2MeanLengthU[int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr][leftGidx:rightGidx+1]])
                        invLengthStr = 'inv: %s%s, %.0fKb' % (invLengthInGenes, 'cs', invLengthInMb * 1000)
                        x_offset = -15
                    else:
                        assert unitOfGenomes == 'gene'
                        invLengthStr = 'inv: %s%s' % (invLengthInGenes, 'cs')
                        x_offset = -7
                    scenarioItems.append(mySvgDrawer.Text(Point(x_offset, interLineFactor * (cptStep-1)), invLengthStr, size=interLineFactor / 3.0))

                    # draw circles around breakpoint reuse
                    leftBreakpointFlanks2 = myIntervals.geneExtremitiesFromIntergeneInChrom(leftGidx, genome1OnlySchr[sChr])
                    rightBreakpointFlanks2 = myIntervals.geneExtremitiesFromIntergeneInChrom(rightGidx + 1, genome1OnlySchr[sChr])
                    isBreakPointReuseL = bra.recordBreakpoint((sChr, leftGidx), genome1OnlySchr)
                    isBreakPointReuseR = bra.recordBreakpoint((sChr, rightGidx + 1), genome1OnlySchr)
                    cpt_nbBreakpointReuse = sum([isBreakPointReuseL, isBreakPointReuseR])
                    nbBreakpointReuse2 += cpt_nbBreakpointReuse
                    if cpt_nbBreakpointReuse > 0:
                        if isBreakPointReuseL:
                            scenarioItems.append(mySvgDrawer.Circle(Point(leftLengthGraph, interLineFactor * (cptStep-1)), radiusBR, (255, 0, 0), opacity=BRopacity))
                            for geneExtr in myIntervals.geneExtremitiesFromIntergeneInChrom(leftGidx, genome1OnlySchr[sChr]):
                                if geneExtr in setRemnantBreakpointFlanks2:
                                    dictBreakpointReusePerFlank[geneExtr] += 1
                        if isBreakPointReuseR:
                            scenarioItems.append(mySvgDrawer.Circle(Point(rightLengthGraph, interLineFactor * (cptStep-1)), radiusBR, (255, 0, 0), opacity=BRopacity))
                            for geneExtr in myIntervals.geneExtremitiesFromIntergeneInChrom(rightGidx + 1, genome1OnlySchr[sChr]):
                                if geneExtr in setRemnantBreakpointFlanks2:
                                    dictBreakpointReusePerFlank[geneExtr] += 1
                        # write the nb of breakpoints reused on the right
                        scenarioItems.append(mySvgDrawer.Text(Point(101, interLineFactor * (cptStep-1)), '%s B.R.' % cpt_nbBreakpointReuse, size=interLineFactor / 3.0))
                    setRemnantBreakpointFlanks2 |= leftBreakpointFlanks2 | rightBreakpointFlanks2
                    ############################################
                    # perform the inversions
                    ############################################
                    genome1OnlySchr = myEvents.performInversion(genome1OnlySchr, (sChr, leftGidx, rightGidx + 1))
                    ############################################
                    bra.recordNewAdjacencyFromNewIntergene(genome1OnlySchr, (sChr, leftGidx))
                    bra.recordNewAdjacencyFromNewIntergene(genome1OnlySchr, (sChr, rightGidx + 1))

                    # change cs names
                    genome1OnlySchrWithOldNames = myLightGenomes.LightGenome(genome1OnlySchr)
                    genome1OnlySchrWithOldNames[sChr] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome1OnlySchr[sChr]]
                    symbolsInGenes = [g.n for g in genome1OnlySchrWithOldNames[sChr]]
                    lengthGenes = [sbId2MeanLengthU[int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr]]
                    (chromosomeItems, factorLength2) = myGenomesDrawer.svgItemsChromosome(genome1OnlySchrWithOldNames[sChr],
                                                                                      lengthGenes=lengthGenes,
                                                                                      halfIntergeneLengths=halfIntergeneLength,
                                                                                      desiredLength=100,
                                                                                      symbolsInGenes=symbolsInGenes,
                                                                                      geneWidth=widthGene,
                                                                                      name2color=name2color)
                    assert abs(factorLength - factorLength2) < 0.0001
                    scenarioItems += mySvgDrawer.translateItems(chromosomeItems, (0, interLineFactor * cptStep))
                    # print >> sys.stderr, [('-' if s == -1 else '') + oldGn2newGn.inv[n] for (n, s) in genome1OnlySchr[sChr]]
            elif stepWasOnPrevLine:
                stepWasOnPrevLine = False
            else:
                stepWasOnPrevLine = False
        # freeze the dictdefault to avoid any future mistake
        dictBreakpointReusePerFlank = dict(dictBreakpointReusePerFlank)
        #assert set(dictBreakpointReusePerFlank.keys()) == setGeneExtremitiesWithBreakpointReuse, '%s, %s' % (str(set(dictBreakpointReusePerFlank.keys()) - setGeneExtremitiesWithBreakpointReuse),  str(setGeneExtremitiesWithBreakpointReuse - set(dictBreakpointReusePerFlank.keys())))
        #assert nbBreakpointReuse1 == nbBreakpointReuse2
        print >> sys.stderr, 'nb breakpoints = %s' % nbBreakpoints
        print >> sys.stderr, 'nb breakpoints reuse = %s' % nbBreakpointReuse2
        print >> sys.stderr, 'nb breakpoints at chromosomes extremities = %s' % nbBreakpointsAtChromExtremity

        print >> sys.stderr, '############# bra statistics ##########'
        bra.printStats()
        # next line can be used to know the remaining adjancencies between two sbs that were not broken
        # print >> sys.stderr, bra.conservedAncestralBlocks

        # scenarioItems = mySvgDrawer.zItems(scenarioItems, mySvgDrawer.Gene)
        scenarioItems = mySvgDrawer.translateItems(scenarioItems, (0, yOrigin))
        sceneItems += scenarioItems

        # scale the mb to the scale of the first genome
        sbCoordsInUChr2 = [(l*lengthFactor, r*lengthFactor) for (l, r) in sbCoordsInUChr2]

        symbolsInGenes = [g.n for g in genome2OnlySchrWithOldNames[sChr]]
        items_genomeWithSbCoords2 = myGenomesDrawer.svgItemsChromosome(genome2OnlySchrWithOldNames[sChr],
                                                                       colorsGenerator=None,
                                                                       genesCoordinates=sbCoordsInUChr2,
                                                                       lengthGenes=None,
                                                                       halfIntergeneLengths=None,
                                                                       desiredLength=None,
                                                                       symbolsInGenes=symbolsInGenes,
                                                                       geneWidth=widthGene,
                                                                       name2color=name2color)


        # draw localisation of breakpoints reused on the mouse chromosome
        # FIXME, for the moment both telomeres have the same lengths
        halfTelomere = sbCoordsInUChr2[0][0] / 2.0
        # 1st telomere
        genome2OnlySchr.computeDictG2Ps()
        def radiusOfNbBreakPointReuse(breakPointReuse, maxAllowedCircleRadius=(interLineFactor/1.5), minCircleRadius=radiusBR):
            maxAllowedCircleSurface = math.pi * pow(maxAllowedCircleRadius, 2)
            maxBreakpointReuse = max(dictBreakpointReusePerFlank.values())
            slope = (float(maxAllowedCircleSurface) - minCircleRadius) / maxBreakpointReuse
            surface = minCircleRadius + slope * breakPointReuse
            radius = math.sqrt(surface / math.pi)
            return radius

        lge = myIntervals.geneExtremityFromGene(genome2OnlySchr[sChr][0], -1)
        if lge in dictBreakpointReusePerFlank:
            loc = halfTelomere
            _radiusBR = radiusOfNbBreakPointReuse(dictBreakpointReusePerFlank[lge])
            items_genomeWithSbCoords2.append(mySvgDrawer.Circle(Point(loc, 0), _radiusBR, (255, 0, 0), opacity=BRopacity))
            old_loc = loc
        for (idxG, oadj) in enumerate(myIntervals.analyseGenomeIntoAdjacencies(genome2OnlySchr, oriented=True, asA=list, fixOrderInAdj=True)):
            (lge, rge) = myIntervals.geneExtremitiesFromAdjacency(oadj)
            if lge in dictBreakpointReusePerFlank or rge in dictBreakpointReusePerFlank:
                lgePositions = genome2OnlySchr.getPositions(lge.n)
                assert len(lgePositions) == 1
                rgePositions = genome2OnlySchr.getPositions(rge.n)
                assert len(rgePositions) == 1
                assert next(iter(lgePositions)).idx < next(iter(rgePositions)).idx
                assert sbCoordsInUChr2[idxG][0] < sbCoordsInUChr2[idxG][1]
                loc_left = sbCoordsInUChr2[idxG][1]
                assert sbCoordsInUChr2[idxG+1][0] < sbCoordsInUChr2[idxG+1][1]
                loc_right = sbCoordsInUChr2[idxG+1][0]
                if not (loc_left < loc_right):
                    print >> sys.stderr, ' due to a small overlap, the positioning of breakpoints reused is not perfect in the last genome, %s < %s' % (loc_left, loc_right)
                loc = (loc_left + loc_right) / 2.0
                nbBreakPointReuse =  max(dictBreakpointReusePerFlank[lge] if lge in dictBreakpointReusePerFlank else 0,
                                         dictBreakpointReusePerFlank[rge] if rge in dictBreakpointReusePerFlank else 0)
                _radiusBR = radiusOfNbBreakPointReuse(nbBreakPointReuse)
                items_genomeWithSbCoords2.append(mySvgDrawer.Circle(Point(loc, 0), _radiusBR, (255, 0, 0), opacity=BRopacity))
                old_loc = loc
        rge = myIntervals.geneExtremityFromGene(genome2OnlySchr[sChr][-1], +1)
        # last telomere
        if rge in dictBreakpointReusePerFlank:
            loc = sbCoordsInUChr2[-1][1] + halfTelomere
            assert old_loc <= loc
            _radiusBR = radiusOfNbBreakPointReuse(dictBreakpointReusePerFlank[rge])
            items_genomeWithSbCoords2.append(mySvgDrawer.Circle(Point(loc, 0), _radiusBR, (255, 0, 0), opacity=BRopacity))

        lengthGenome2 = max([l.end for l in items_genomeWithSbCoords2 if isinstance(l, mySvgDrawer.Line)]).x
        # translate in the x_axis to center
        txGenome2 = (100 - lengthGenome2) / 2.0
        items_genomeWithSbCoords2 = mySvgDrawer.translateItems(items_genomeWithSbCoords2, (txGenome2, yOrigin + (cptStep + 1) * interLineFactor))

        sceneItems += items_genomeWithSbCoords2
        (origin, wwidth, hheight) = mySvgDrawer.boundingBoxItems(sceneItems)
        origin = origin + Point(-5, 0)
        wwidth += 5
        print >> sys.stderr, str((origin, wwidth, hheight))
        # (interLineFactor * (cptStep + 3))
        scene = mySvgDrawer.Scene(name='chromosome', origin=origin, width=wwidth, height=hheight)
        for item in sceneItems:
            scene.add(item)
        scene.write_svg(filename=arguments['outPath'] + '/svgScenario.svg')

        print >> seevolutionInputFile, '</speciation>\n</organism>\n</genomemutations>'
        seevolutionInputFile.close()
    else:
        raise ValueError('The path to the executable file of grimm is incorrect')

if __name__ == '__main__':
    main()