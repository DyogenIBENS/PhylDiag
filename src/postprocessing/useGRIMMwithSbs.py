#!/usr/bin/python
# -*- coding: utf-8 -*-

# This script wraps the GRIMM v2.01 software
# http://grimm.ucsd.edu/GRIMM/

# The aim is to use PhylDiag's conserved segments between human and mouse in the X chromosome
# and use these conserved segments to retrieve reversal scenarios infered from GRIMM

import collections
import itertools
from numpy import random
import sys
import os
import subprocess
from utils import myTools, myDiags, myLightGenomes

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

def rewriteGenomesIntoSbs(sbsInPairComp, checkNoOverlap=True):
    assert isinstance(sbsInPairComp, myTools.Dict2d)

    # give ids to each sbs
    sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        sbsInPairCompWithIds.addToLocation((c1, c2), sb)

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
                                  [('outPath', str, './res'), ('pathToGRIMM', str, '/home/jlucas/Libs/GRIMM-2.01/grimm')],
                                  __doc__)
    genome1 = myLightGenomes.LightGenome(arguments['genome1'], withDict=True)
    genome2 = myLightGenomes.LightGenome(arguments['genome2'], withDict=True)
    sbsInPairComp = myDiags.parseSbsFile(arguments['syntenyBlocks'],
                                         genome1=genome1,
                                         genome2=genome2)

    (genome1InSbs, genome2InSbs) = rewriteGenomesIntoSbs(sbsInPairComp)

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
    genome1WithOnlyX['X'] = genome1InSbs['X']
    genome2WithOnlyX['X'] = genome2InSbs['X']
    setGeneNames = genome1WithOnlyX.getGeneNames() & genome2WithOnlyX.getGeneNames()
    genome1WithOnlyX.removeGenes(genome1WithOnlyX.getGeneNames() - setGeneNames)
    genome1WithOnlyX.removeGenes(genome2WithOnlyX.getGeneNames() - setGeneNames)
    dictOldGeneNameToNewGeneName = {}
    cptGeneName = 1
    for gn in setGeneNames:
        dictOldGeneNameToNewGeneName[gn] = str(cptGeneName)
        cptGeneName += 1
    for i, g in enumerate(genome1WithOnlyX['X']):
        genome1WithOnlyX['X'][i] = myLightGenomes.OGene(dictOldGeneNameToNewGeneName[g.n], g.s)
    for i, g in enumerate(genome2WithOnlyX['X']):
        genome2WithOnlyX['X'][i] = myLightGenomes.OGene(dictOldGeneNameToNewGeneName[g.n], g.s)

    pathToSyntenyBlocksForGrimm = arguments['outPath'] + '/syntenyBlocksForGRIMM.txt'
    with open(pathToSyntenyBlocksForGrimm, 'w') as f:
        printGenomeInto_GRIMM_format(genome1WithOnlyX, stream=f)
        printGenomeInto_GRIMM_format(genome2WithOnlyX, stream=f)

    # wrap the grimm executable
    pathToGrimm = arguments['pathToGRIMM']
    if os.path.isfile(pathToGrimm):
        bashComand = pathToGrimm + ' -f ' + pathToSyntenyBlocksForGrimm
        subProcStdin = None  #  subprocess.PIPE, if you want it
        subProcStderr = None  #  subprocess.PIPE, if you want it
        proc = subprocess.Popen(bashComand, shell=True, stdin=subProcStdin, stdout=subprocess.PIPE, stderr=subProcStderr)
        # returns stdout and stderr
        stdout, subProcStderr = proc.communicate(subProcStdin)
        # proc.returncode, might be interesting too
        print >> sys.stdout, stdout
    else:
        raise ValueError('The path to the executable file of grimm is incorrect')

if __name__ == '__main__':
    main()