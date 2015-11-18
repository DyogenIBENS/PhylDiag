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
                assert min(sb1.maxOnG(1)) < max(sb2.minOnG(1)), "%s < %s" % (sb1.maxOnG(1), sb2.minOnG(1))
            for ((id1, sb1, s1), (id2, sb2, s2)) in myTools.myIterator.slidingTuple(listOfSortedSbsInC2[c2]):
                assert min(sb1.maxOnG(2)) < max(sb2.minOnG(2)), "%s < %s" % (sb1.maxOnG(2), sb2.minOnG(2))

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
                                   ('renameSbsByOrder', bool, True)],
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
    genome1WithOnlyX = myLightGenomes.LightGenome()
    genome2WithOnlyX = myLightGenomes.LightGenome()
    genome1WithOnlyX.name = genome1.name
    genome2WithOnlyX.name = genome2.name
    genome1WithOnlyX['X'] = genome1InSbs['X']
    genome2WithOnlyX['X'] = genome2InSbs['X']
    setGeneNames = genome1WithOnlyX.getGeneNames() & genome2WithOnlyX.getGeneNames()
    genome1WithOnlyX.removeGenes(genome1WithOnlyX.getGeneNames() - setGeneNames)
    genome1WithOnlyX.removeGenes(genome2WithOnlyX.getGeneNames() - setGeneNames)

    # renaming should be done because of internal limitations of GRIMM (that do not want gene names to be superior than
    # the total number of genes.
    oldGn2newGn = bidict.bidict()
    if max([int(gn) for (gn, gs) in genome1WithOnlyX['X']]) > len(genome1WithOnlyX['X']):
        cptGeneName = 1
        for gn in setGeneNames:
            oldGn2newGn[gn] = str(cptGeneName)
            cptGeneName += 1
        for i, g in enumerate(genome1WithOnlyX['X']):
            genome1WithOnlyX['X'][i] = myLightGenomes.OGene(oldGn2newGn[g.n], g.s)
        for i, g in enumerate(genome2WithOnlyX['X']):
            genome2WithOnlyX['X'][i] = myLightGenomes.OGene(oldGn2newGn[g.n], g.s)
    else:
        for gn in setGeneNames:
            oldGn2newGn[gn] = gn

    pathToSyntenyBlocksForGrimm = arguments['outPath'] + '/syntenyBlocksForGRIMM.txt'
    with open(pathToSyntenyBlocksForGrimm, 'w') as f:
        printGenomeInto_GRIMM_format(genome1WithOnlyX, stream=f)
        printGenomeInto_GRIMM_format(genome2WithOnlyX, stream=f)
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
        with open('res/grimmOutputFile.txt', 'w') as grimmOutputFile:
            print >> grimmOutputFile, grimmStdOut

        approx_equal = lambda a, b, t: abs(a - b) < t
        stepWasOnPrevLine = False
        # print >> sys.stderr, [oldGn2newGn.inv[gn] for (gn, s) in genome1WithOnlyX['X']]
        # print >> sys.stderr, sbId2Length.keys()
        totalLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1WithOnlyX['X']])
        import libs.myEvents as myEvents
        for (n, s) in genome1WithOnlyX['X']:
            assert s == +1 or s == -1

        seevolutionInputFile = open('res/grimmOutputFile.txt', 'w')
        # input xml file for seevolution
        print >> seevolutionInputFile, '<?xml version="1.0" encoding="ISO-8859-1"?>'
        print >> seevolutionInputFile, '<genomemutations id="1">'
        print >> seevolutionInputFile, '<organism id="Organism">'
        # the small +1 is here to avoid bound error in seevolution
        # only the circular version works
        print >> seevolutionInputFile, '<chromosome id="%s" circular="true" length="%s"/>' % ('X', int(round(totalLength*1000, 0)))
        print >> seevolutionInputFile, '<speciation name="ancestor">'

        cptStep = 0
        scenarioItems = []
        halfIntergeneLength = 1.5
        interLineFactor = 6.0
        import libs.myBreakpointsAnalyser
        bra = libs.myBreakpointsAnalyser.BreakpointsAnalyser(genome1WithOnlyX)
        nbBreakpoints = 0
        nbBreakpointReuse = 0
        nbBreakpointsAtChromExtremity = 0
        setBreakpointFlanks = set([])
        for line in grimmStdOut.split('\n'):
            print >> sys.stderr, line
            if 'Step' in line:
                stepWasOnPrevLine = True
                if '(Source)' in line:
                    genome1WithOnlyXWithOldNames = myLightGenomes.LightGenome(genome1WithOnlyX)
                    genome1WithOnlyXWithOldNames['X'] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome1WithOnlyX['X']]
                    scenarioItems += myGenomesDrawer.drawLightGenome(genome1WithOnlyXWithOldNames,
                                                                  lengthGenes={'X': [sbId2Length[int(gn)] for (gn, s) in genome1WithOnlyXWithOldNames['X']]},
                                                                  halfIntergeneLength=halfIntergeneLength)['X']
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
                    # print >> sys.stderr, '(leftGn, rightGn)= (%s, %s)' % (leftGn, rightGn)
                    assert genome1WithOnlyX['X'][leftGidx].n == leftGn, '%s == %s' % (genome1WithOnlyX['X'][leftGidx].n, leftGn)
                    assert genome1WithOnlyX['X'][rightGidx].n == rightGn, '%s == %s' % (genome1WithOnlyX['X'][rightGidx].n, rightGn)

                    leftLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1WithOnlyX['X'][:leftGidx]])
                    invLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1WithOnlyX['X'][leftGidx:rightGidx+1]])
                    rightLength = sum([sbId2Length[int(oldGn2newGn.inv[gn])] for (gn, s) in genome1WithOnlyX['X'][rightGidx+1:]])
                    assert approx_equal(leftLength + invLength + rightLength, totalLength, 0.01)

                    # because of seevolution
                    if leftGidx == rightGidx:
                        invLength = 0
                    print >> seevolutionInputFile, "<inversion left=\"%s\" right=\"%s\" chromosome=\"%s\"/>" % (int(round(leftLength*1000, 0)),
                                                                                                      int(round(leftLength*1000, 0)) + int(round(invLength*1000,0)),
                                                                                                      'X')
                    cptStep += 1
                    # analyse breakpoints before performing them
                    nbBreakpoints += 2
                    if leftGidx == 0:
                        nbBreakpointsAtChromExtremity += 1
                    else:
                        assert leftGidx < len(genome1WithOnlyX['X'])
                        leftBreakpointFlanks = myIntervals.geneExtremitiesFromIntergeneInChrom(leftGidx, genome1WithOnlyX['X'])
                        print >> sys.stderr, 'left breakpoint = %s' % leftBreakpointFlanks
                        assert len(leftBreakpointFlanks) == 2
                        if len(leftBreakpointFlanks & setBreakpointFlanks) > 0:
                            nbBreakpointReuse += 1
                    if rightGidx == len(genome1WithOnlyX['X']) - 1:
                        nbBreakpointsAtChromExtremity += 1
                    else:
                        assert rightGidx > 0
                        rightBreakpointFlanks = myIntervals.geneExtremitiesFromIntergeneInChrom(rightGidx + 1, genome1WithOnlyX['X'])
                        print >> sys.stderr, 'right breakpoint = %s' % rightBreakpointFlanks
                        assert len(rightBreakpointFlanks) == 2
                        if len(rightBreakpointFlanks & setBreakpointFlanks) > 0:
                            nbBreakpointReuse += 1
                    setBreakpointFlanks |= leftBreakpointFlanks | rightBreakpointFlanks

                    bra.recordBreakpoint(('X', leftGidx), genome1WithOnlyX)
                    bra.recordBreakpoint(('X', rightGidx + 1), genome1WithOnlyX)
                    # perform the inversions
                    genome1WithOnlyX = myEvents.performInversion(genome1WithOnlyX, ('X', leftGidx, rightGidx + 1))
                    bra.recordNewAdjacencyFromNewIntergene(genome1WithOnlyX, ('X', leftGidx))
                    bra.recordNewAdjacencyFromNewIntergene(genome1WithOnlyX, ('X', rightGidx + 1))

                    # change cs names
                    genome1WithOnlyXWithOldNames = myLightGenomes.LightGenome(genome1WithOnlyX)
                    genome1WithOnlyXWithOldNames['X'] = [myLightGenomes.OGene(oldGn2newGn.inv[gn], s) for (gn, s) in genome1WithOnlyX['X']]
                    genomeItems = myGenomesDrawer.drawLightGenome(genome1WithOnlyXWithOldNames,
                                                                  lengthGenes={'X': [sbId2Length[int(gn)] for (gn, s) in genome1WithOnlyXWithOldNames['X']]},
                                                                  halfIntergeneLength=halfIntergeneLength)
                    scenarioItems += mySvgDrawer.tanslateItems(genomeItems['X'], (0, interLineFactor * cptStep))
                    # 0.1 is the halfintergene length
                    leftLengthGraph = sum([(sbId2Length[int(oldGn2newGn.inv[gn])] + 2 * halfIntergeneLength) for (gn, s) in genome1WithOnlyX['X'][:leftGidx]])
                    invLengthGraph = sum([(sbId2Length[int(oldGn2newGn.inv[gn])] + 2 * halfIntergeneLength) for (gn, s) in genome1WithOnlyX['X'][leftGidx:rightGidx+1]])
                    rightLengthGraph = leftLengthGraph + invLengthGraph
                    leftLengthGraph -= halfIntergeneLength / 2.0
                    rightLengthGraph -= halfIntergeneLength / 2.0
                    # add lines to visualise the inversion
                    scenarioItems.append(mySvgDrawer.Line(Point(leftLengthGraph, interLineFactor * (cptStep-1)), Point(rightLengthGraph, interLineFactor * cptStep), width=0.2))
                    scenarioItems.append(mySvgDrawer.Line(Point(rightLengthGraph, interLineFactor * (cptStep-1)), Point(leftLengthGraph, interLineFactor * cptStep), width=0.2))
                    # print >> sys.stderr, [('-' if s == -1 else '') + oldGn2newGn.inv[n] for (n, s) in genome1WithOnlyX['X']]
            elif stepWasOnPrevLine:
                stepWasOnPrevLine = False
            else:
                stepWasOnPrevLine = False

        # TODO: solve this, the adjacency 1 - 20 is never broken !!!!!!!

        print >> sys.stderr, 'nb breakpoints = %s' % nbBreakpoints
        print >> sys.stderr, 'nb breakpoints reuse = %s' % nbBreakpointReuse
        print >> sys.stderr, 'nb breakpoints at chromosomes extremities = %s' % nbBreakpointsAtChromExtremity

        print >> sys.stderr, '############# bra statistics ##########'
        bra.printStats()

        scenarioItems = mySvgDrawer.zItems(scenarioItems, mySvgDrawer.Gene)
        scene = mySvgDrawer.Scene(name='chromosome', width=100, height=100)
        for item in scenarioItems:
            scene.add(item)
        scene.write_svg(filename='svgScenario.svg')

        print >> seevolutionInputFile, '</speciation>\n</organism>\n</genomemutations>'
        seevolutionInputFile.close()
    else:
        raise ValueError('The path to the executable file of grimm is incorrect')

if __name__ == '__main__':
    main()