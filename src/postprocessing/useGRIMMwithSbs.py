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
from utils import myTools, myDiags, myLightGenomes, myGenomes, myGenomesDrawer, mySvgDrawer, myIntervals
import bidict
from utils. mySvgDrawer import Point as Point
import libs.myEvents as myEvents
import libs.myBreakpointsAnalyser

GENOMENAME = 0

def projectSbsOnGenome(listOfSbsInCompC1C2, cx):
    assert cx in [1, 2]
    listSbsInCx = list(listOfSbsInCompC1C2)
    # give orientations to each sb
    if cx == 1:
        # orientations of blocks are all set to +1 in the reference chromosome
        listSbsInCx = [(id, sb, +1) for (id, sb, dt) in listOfSbsInCompC1C2]
    else:
        assert cx == 2
        for (i, (id, sb, dt)) in enumerate(listSbsInCx):
            if dt == '/':
                s = +1
            elif dt == '\\':
                s = -1
            else:
                assert dt == None
                s = None
            listSbsInCx[i] = (id, sb, s)
    return listSbsInCx

def rewriteGenomesIntoSbs(sbsInPairComp, checkNoOverlap=True, inputWithIds=True):

    if inputWithIds == False:
        # give ids to each sbs
        assert isinstance(sbsInPairComp, myTools.Dict2d)
        sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
        for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
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
        listOfSbsInCompC1C2 = []
        for idx, sb in enumerate(sbsInPairCompWithIds[c1][c2]):
            # find diagType
            dt = sb.dt
            listOfSbsInCompC1C2.append((sbsInPairCompWithIds.location2id[c1][c2][idx], sb, dt))
        listOfSortedSbsInC1[c1] += projectSbsOnGenome(listOfSbsInCompC1C2, 1)
        listOfSortedSbsInC2[c2] += projectSbsOnGenome(listOfSbsInCompC1C2, 2)

    for (c1, c2) in [(c1, c2) for c1 in sbsInPairCompWithIds for c2 in sbsInPairCompWithIds[c1]]:
        # sort sbs by starting position on each genome
        listOfSortedSbsInC1[c1].sort(key=lambda x: x[1].minOnG(1))
        listOfSortedSbsInC2[c2].sort(key=lambda x: x[1].minOnG(2))

        if checkNoOverlap:
            for ((id1, sb1, s1), (id2, sb2, s2)) in myTools.myIterator.slidingTuple(listOfSortedSbsInC1[c1]):
                if sb2.minOnG(1) < sb1.maxOnG(1):
                    print >> sys.stderr, 'Warning, overlap of sbs:' + "%s < %s" % (sb1.maxOnG(1), sb2.minOnG(1))
                # assert sb1.maxOnG(1) < sb2.minOnG(1), "%s < %s" % (sb1.maxOnG(1), sb2.minOnG(1))
            for ((id1, sb1, s1), (id2, sb2, s2)) in myTools.myIterator.slidingTuple(listOfSortedSbsInC2[c2]):
                #assert sb1.maxOnG(2) < sb2.minOnG(2), "%s < %s" % (sb1.maxOnG(2), sb2.minOnG(2))
                if sb2.minOnG(2) < sb1.maxOnG(2):
                    print >> sys.stderr, 'Warning, overlap of sbs:' + "%s < %s" % (sb1.maxOnG(2), sb2.minOnG(2))

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
                                  [('outPath', str, './res'),
                                   ('pathToGRIMM', str, '/home/jlucas/Libs/GRIMM-2.01/grimm'),
                                   ('renameSbsByOrder', bool, True),
                                   ('selectedChrom', str, 'X')],
                                  __doc__)
    genome1 = myLightGenomes.LightGenome(arguments['genome1'], withDict=True)
    genome2 = myLightGenomes.LightGenome(arguments['genome2'], withDict=True)
    sbsInPairCompWithIds = myDiags.parseSbsFile(arguments['syntenyBlocks'],
                                         genome1=genome1,
                                         genome2=genome2,
                                         withIds=True)
    # rename : id=1 (longest sb), .....
    if arguments['renameSbsByOrder']:
        new_sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
        newId = 1
        for (id, (k1, k2), sb) in sorted([(id, (k1, k2), sb) for (id, (k1, k2), sb) in sbsInPairCompWithIds.iterByOrderedIds()], key=lambda x: len(x[2].la), reverse=True):
            new_sbsInPairCompWithIds.addToLocationWithId((k1, k2), sb, newId)
            newId += 1
        sbsInPairCompWithIds = new_sbsInPairCompWithIds
        with open(arguments['outPath'] + '/syntenyBlocksNewNames.txt', 'w') as f:
            myDiags.printSbsFile(sbsInPairCompWithIds, genome1, genome2, stream=f)

    # get the lengths of sbs
    genome1H = myGenomes.Genome(arguments['genome1'])
    genome2H = myGenomes.Genome(arguments['genome2'])
    sbId2Length = myDiags.getSbsLengths(genome1H, genome2H, sbsInPairCompWithIds, lengthUnit='Mb')

    # FIXME
    sbIdCoordsInGenome1 = collections.defaultdict(lambda: collections.defaultdict(tuple))
    sbIdCoordsInGenome2 = collections.defaultdict(lambda: collections.defaultdict(tuple))
    def findSbCoordOnGenome(sb, c, genomeH, rankGenome):
        assert isinstance(sb, myDiags.SyntenyBlock)
        (minIdx, maxIdx) = (sb.minOnG(rankGenome), sb.maxOnG(rankGenome))
        assert isinstance(minIdx, int)
        leftGene = genomeH.lstGenes[c][minIdx]
        assert isinstance(leftGene, myGenomes.Gene)
        rightGene = genomeH.lstGenes[c][maxIdx]
        return (leftGene.beginning / 1000000.0, rightGene.end / 1000000.0)

    for (id, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        sbIdCoordsInGenome1[c1][id] = findSbCoordOnGenome(sb, myGenomes.commonChrName(c1), genome1H, 1)
        sbIdCoordsInGenome2[c2][id] = findSbCoordOnGenome(sb, myGenomes.commonChrName(c2), genome2H, 2)

    (genome1InSbs, genome2InSbs) = rewriteGenomesIntoSbs(sbsInPairCompWithIds)
    genome1InSbs.name = genome1.name
    genome2InSbs.name = genome2.name
    genome1InSbs.sort(byName=True)
    genome2InSbs.sort(byName=True)
    # print these two representation of the genomes
    pathToSyntenyBlocksForGrimm = arguments['outPath'] + '/syntenyBlocksForGRIMM_allChromosomes.txt'
    with open(pathToSyntenyBlocksForGrimm, 'w') as f:
        printGenomeInto_GRIMM_format(genome1InSbs, stream=f)
        printGenomeInto_GRIMM_format(genome2InSbs, stream=f)

    # reduce genomes to only chromosomes X
    genome1OnlySchr = myLightGenomes.LightGenome()
    genome2OnlySchr = myLightGenomes.LightGenome()
    genome1OnlySchr.name = genome1.name
    genome2OnlySchr.name = genome2.name
    sChr = arguments['selectedChrom']
    genome1OnlySchr[sChr] = genome1InSbs[sChr]
    genome2OnlySchr[sChr] = genome2InSbs[sChr]
    setGeneNames = genome1OnlySchr.getGeneNames() & genome2OnlySchr.getGeneNames()
    genome1OnlySchr.removeGenes(genome1OnlySchr.getGeneNames() - setGeneNames)
    genome2OnlySchr.removeGenes(genome2OnlySchr.getGeneNames() - setGeneNames)

    # renaming should be done because of internal limitations of GRIMM (that do not want gene names to be superior than
    # the total number of genes.
    oldGn2newGn = bidict.bidict()
    if max([int(gn) for (gn, gs) in genome1OnlySchr[sChr]]) > len(genome1OnlySchr[sChr]):
        cptGeneName = 1
        for gn in setGeneNames:
            oldGn2newGn[gn] = str(cptGeneName)
            cptGeneName += 1
        for i, g in enumerate(genome1OnlySchr[sChr]):
            genome1OnlySchr[sChr][i] = myLightGenomes.OGene(oldGn2newGn[g.n], g.s)
        for i, g in enumerate(genome2OnlySchr[sChr]):
            genome2OnlySchr[sChr][i] = myLightGenomes.OGene(oldGn2newGn[g.n], g.s)
    else:
        for gn in setGeneNames:
            oldGn2newGn[gn] = gn

    pathToSyntenyBlocksForGrimm = arguments['outPath'] + '/syntenyBlocksForGRIMM.txt'
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
        totalLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr]])
        for (n, s) in genome1OnlySchr[sChr]:
            assert s == +1 or s == -1

        seevolutionInputFile = open(arguments['outPath'] + '/grimmOutputFile.txt', 'w')
        # input xml file for seevolution
        print >> seevolutionInputFile, '<?xml version="1.0" encoding="ISO-8859-1"?>'
        print >> seevolutionInputFile, '<genomemutations id="1">'
        print >> seevolutionInputFile, '<organism id="Organism">'
        # the small +1 is here to avoid bound error in seevolution
        # only the circular version works
        print >> seevolutionInputFile, '<chromosome id="%s" circular="true" length="%s"/>' % (sChr, int(round(totalLength*1000, 0)))
        print >> seevolutionInputFile, '<speciation name="ancestor">'

        # Draw the human genome with lengths in the human
        genome1OnlySchrWithOldNames = myLightGenomes.LightGenome(genome1OnlySchr)
        genome1OnlySchrWithOldNames[sChr] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome1OnlySchr[sChr]]

        genome2OnlySchrWithOldNames = myLightGenomes.LightGenome(genome2OnlySchr)
        genome2OnlySchrWithOldNames[sChr] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome2OnlySchr[sChr]]

        # print >> sys.stderr, 'sbIdCoordsInGenome1[\'X\']=', sbIdCoordsInGenome1[sChr]
        sbCoordsX1 = {sChr: [sbIdCoordsInGenome1[sChr][int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr]]}
        # print >> sys.stderr, 'sbCoordsX1=', sbCoordsX1
        sbCoordsX2 = {sChr: [sbIdCoordsInGenome2[sChr][int(gn)] for (gn, s) in genome2OnlySchrWithOldNames[sChr]]}

        (items_genomeWithSbCoords1, lengthFactorBp) = myGenomesDrawer.svgItemsLightGenome(genome1OnlySchrWithOldNames,
                                                                                          colorsGenerator=None,
                                                                                          genesCoordinates=sbCoordsX1,
                                                                                          lengthGenes=None,
                                                                                          desiredLength=100)
        items_genomeWithSbCoords1 = items_genomeWithSbCoords1[sChr]

        cptStep = 1
        scenarioItems = []
        halfIntergeneLength = 1.5
        radiusBR = 0.5
        BRopacity = 0.5
        interLineFactor = 3.0
        yOrigin = 1 * interLineFactor
        sceneItems = mySvgDrawer.translateItems(items_genomeWithSbCoords1, (0, yOrigin))

        bra = libs.myBreakpointsAnalyser.BreakpointsAnalyser(genome1OnlySchr)
        setGeneExtrBreakReuse = collections.defaultdict(int)
        nbBreakpoints = 0
        nbBreakpointReuse = 0
        nbBreakpointsAtChromExtremity = 0
        setBreakpointFlanks = set([])
        for line in grimmStdOut.split('\n'):
            print >> sys.stderr, line
            if 'Step' in line:
                stepWasOnPrevLine = True
                if '(Source)' in line:
                    scenarioItems = myGenomesDrawer.svgItemsLightGenome(genome1OnlySchrWithOldNames,
                                                                        lengthGenes={sChr: [sbId2Length[int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr]]},
                                                                        halfIntergeneLengths=halfIntergeneLength,
                                                                        desiredLength=100)[0][sChr]
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

                    leftLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr][:leftGidx]])
                    invLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr][leftGidx:rightGidx+1]])
                    rightLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1OnlySchr[sChr][rightGidx+1:]])
                    assert approx_equal(leftLength + invLength + rightLength, totalLength, 0.01)

                    # because of seevolution
                    if leftGidx == rightGidx:
                        invLength = 0
                    print >> seevolutionInputFile, "<inversion left=\"%s\" right=\"%s\" chromosome=\"%s\"/>" % (int(round(leftLength*1000, 0)),
                                                                                                      int(round(leftLength*1000, 0)) + int(round(invLength*1000,0)),
                                                                                                      sChr)
                    cptStep += 1
                    # analyse breakpoints before performing them
                    nbBreakpoints += 2
                    rightBreakpointFlanks = set()
                    leftBreakpointFlanks = set()
                    if leftGidx == 0:
                        nbBreakpointsAtChromExtremity += 1
                    else:
                        assert leftGidx < len(genome1OnlySchr[sChr])
                        leftBreakpointFlanks = myIntervals.geneExtremitiesFromIntergeneInChrom(leftGidx, genome1OnlySchr[sChr])
                        assert len(leftBreakpointFlanks) == 2
                        if len(leftBreakpointFlanks & setBreakpointFlanks) > 0:
                            nbBreakpointReuse += 1
                    if rightGidx == len(genome1OnlySchr[sChr]) - 1:
                        nbBreakpointsAtChromExtremity += 1
                    else:
                        print >> sys.stderr, rightGidx
                        assert rightGidx >= 0
                        rightBreakpointFlanks = myIntervals.geneExtremitiesFromIntergeneInChrom(rightGidx + 1, genome1OnlySchr[sChr])
                        assert len(rightBreakpointFlanks) == 2
                        if len(rightBreakpointFlanks & setBreakpointFlanks) > 0:
                            nbBreakpointReuse += 1
                    setBreakpointFlanks |= leftBreakpointFlanks | rightBreakpointFlanks

                    isBreakPointReuseL = bra.recordBreakpoint((sChr, leftGidx), genome1OnlySchr)
                    isBreakPointReuseR = bra.recordBreakpoint((sChr, rightGidx + 1), genome1OnlySchr)
                    # perform the inversions
                    genome1OnlySchr = myEvents.performInversion(genome1OnlySchr, (sChr, leftGidx, rightGidx + 1))
                    bra.recordNewAdjacencyFromNewIntergene(genome1OnlySchr, (sChr, leftGidx))
                    bra.recordNewAdjacencyFromNewIntergene(genome1OnlySchr, (sChr, rightGidx + 1))

                    # change cs names
                    genome1OnlySchrWithOldNames = myLightGenomes.LightGenome(genome1OnlySchr)
                    genome1OnlySchrWithOldNames[sChr] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome1OnlySchr[sChr]]
                    (genomeItems, factorLength) = myGenomesDrawer.svgItemsLightGenome(genome1OnlySchrWithOldNames,
                                                                                      lengthGenes={sChr: [sbId2Length[int(gn)] for (gn, s) in genome1OnlySchrWithOldNames[sChr]]},
                                                                                      halfIntergeneLengths=halfIntergeneLength,
                                                                                      desiredLength=100)
                    nbBreakpointReuse = sum([isBreakPointReuseL, isBreakPointReuseR])
                    scenarioItems += mySvgDrawer.translateItems(genomeItems[sChr], (0, interLineFactor * cptStep))
                    # 0.1 is the halfintergene length
                    leftLengthGraph = factorLength * sum([(sbId2Length[int(gn)] + 2 * halfIntergeneLength) for (gn, s) in genome1OnlySchrWithOldNames[sChr][:leftGidx]])
                    invLengthGraph = factorLength * sum([(sbId2Length[int(gn)] + 2 * halfIntergeneLength) for (gn, s) in genome1OnlySchrWithOldNames[sChr][leftGidx:rightGidx+1]])
                    rightLengthGraph = leftLengthGraph + invLengthGraph
                    # leftLengthGraph -= halfIntergeneLength / 2.0
                    # rightLengthGraph -= halfIntergeneLength / 2.0
                    # add lines to visualise the inversion
                    scenarioItems.append(mySvgDrawer.Line(Point(leftLengthGraph, interLineFactor * (cptStep-1)), Point(rightLengthGraph, interLineFactor * cptStep), width=0.2))
                    scenarioItems.append(mySvgDrawer.Line(Point(rightLengthGraph, interLineFactor * (cptStep-1)), Point(leftLengthGraph, interLineFactor * cptStep), width=0.2))
                    if nbBreakpointReuse > 0:
                        if isBreakPointReuseL:
                            scenarioItems.append(mySvgDrawer.Circle(Point(leftLengthGraph, interLineFactor * (cptStep-1)), radiusBR, (255, 0, 0), opacity=BRopacity))
                            for geneExtr in myIntervals.geneExtremitiesFromIntergeneInChrom(leftGidx, genome1OnlySchr[sChr]):
                                setGeneExtrBreakReuse[geneExtr] += 1
                        if isBreakPointReuseR:
                            scenarioItems.append(mySvgDrawer.Circle(Point(rightLengthGraph, interLineFactor * (cptStep-1)), radiusBR, (255, 0, 0), opacity=BRopacity))
                            for geneExtr in myIntervals.geneExtremitiesFromIntergeneInChrom(rightGidx + 1, genome1OnlySchr[sChr]):
                                setGeneExtrBreakReuse[geneExtr] += 1
                        scenarioItems.append(mySvgDrawer.Text(Point(101, interLineFactor * (cptStep-1)), '%s B.R.' % nbBreakpointReuse, size=interLineFactor / 3.0))
                    # print >> sys.stderr, [('-' if s == -1 else '') + oldGn2newGn.inv[n] for (n, s) in genome1OnlySchr[sChr]]
            elif stepWasOnPrevLine:
                stepWasOnPrevLine = False
            else:
                stepWasOnPrevLine = False

        print >> sys.stderr, 'nb breakpoints = %s' % nbBreakpoints
        print >> sys.stderr, 'nb breakpoints reuse = %s' % nbBreakpointReuse
        print >> sys.stderr, 'nb breakpoints at chromosomes extremities = %s' % nbBreakpointsAtChromExtremity

        print >> sys.stderr, '############# bra statistics ##########'
        bra.printStats()
        # next line can be used to know the remaining adjancencies between two sbs that were not broken
        # print >> sys.stderr, bra.conservedAncestralBlocks

        # scenarioItems = mySvgDrawer.zItems(scenarioItems, mySvgDrawer.Gene)
        scenarioItems = mySvgDrawer.translateItems(scenarioItems, (0, yOrigin))
        sceneItems += scenarioItems

        sbCoordsX2[sChr] = [(l*lengthFactorBp, r*lengthFactorBp) for (l, r) in sbCoordsX2[sChr]]
        items_genomeWithSbCoords2 = myGenomesDrawer.svgItemsLightGenome(genome2OnlySchrWithOldNames,
                                                                        colorsGenerator=None,
                                                                        genesCoordinates=sbCoordsX2,
                                                                        lengthGenes=None,
                                                                        desiredLength=None)[sChr]
        # draw localisation of breakpoints reused on the mouse chromosome
        # FIXME, for the moment both telomeres have the same length
        halfTelomere = sbCoordsX2[sChr][0][0] / 2.0
        # 1st telomere
        lge = myIntervals.geneExtremityFromGene(genome2OnlySchr[sChr][0], -1)
        if lge in setGeneExtrBreakReuse:
            loc = halfTelomere
            items_genomeWithSbCoords2.append(mySvgDrawer.Circle(Point(loc, 0), radiusBR, (255, 0, 0), opacity=BRopacity))
        for (idxG, oadj) in enumerate(myIntervals.analyseGenomeIntoAdjacencies(genome2OnlySchr, oriented=True, asA=list)):
            (lge, rge) = myIntervals.geneExtremitiesFromAdjacency(oadj)
            if lge in setGeneExtrBreakReuse or rge in setGeneExtrBreakReuse:
                (loc_left, loc_right) = (sbCoordsX2[sChr][idxG][1], sbCoordsX2[sChr][idxG+1][0])
                loc = (loc_left + loc_right) / 2.0
                items_genomeWithSbCoords2.append(mySvgDrawer.Circle(Point(loc, 0), radiusBR, (255, 0, 0), opacity=BRopacity))
        lge = myIntervals.geneExtremityFromGene(genome2OnlySchr[sChr][-1], +1)
        # last telomere
        if lge in setGeneExtrBreakReuse:
            loc = sbCoordsX2[sChr][-1][1] + halfTelomere
            items_genomeWithSbCoords2.append(mySvgDrawer.Circle(Point(loc, 0), radiusBR, (255, 0, 0), opacity=BRopacity))

        lengthGenome2 = max([l.end for l in items_genomeWithSbCoords2 if isinstance(l, mySvgDrawer.Line)]).x
        # translate in the x_axis to center
        txGenome2 = (100 - lengthGenome2) / 2.0
        items_genomeWithSbCoords2 = mySvgDrawer.translateItems(items_genomeWithSbCoords2, (txGenome2, yOrigin + (cptStep + 1) * interLineFactor))

        sceneItems += items_genomeWithSbCoords2

        scene = mySvgDrawer.Scene(name='chromosome', width=110, height=(interLineFactor * (cptStep + 3)))
        for item in sceneItems:
            scene.add(item)
        scene.write_svg(filename=arguments['outPath'] + '/svgScenario.svg')

        print >> seevolutionInputFile, '</speciation>\n</organism>\n</genomemutations>'
        seevolutionInputFile.close()
    else:
        raise ValueError('The path to the executable file of grimm is incorrect')

if __name__ == '__main__':
    main()