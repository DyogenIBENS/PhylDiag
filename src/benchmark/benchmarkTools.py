#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import collections
import utils.myMapping as myMapping
import utils.myIntervals as myIntervals
import utils.myLightGenomes as myLightGenomes
from utils import myIntervals

__doc__ = """Tools for benchmarking PhylDiag, i-ADHoRe 3.0 and Cyntenator based on simulated conserved segments (c.s.)"""

def analyseMistakeOfPhylDiag(geneNameInvolvedInAMistake, sbsInPairComp, simulatedSbsGenome,
                             genome1, genome2, mistake='Fp'):
    gn = geneNameInvolvedInAMistake
    assert mistake in {'Fp', 'Fn'}
    assert isinstance(simulatedSbsGenome, myLightGenomes.LightGenome)
    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    assert genome1.withDict is True
    assert genome2.withDict is True
    assert simulatedSbsGenome.withDict is True
    positionsOnGenomes = []
    if mistake == 'Fp':
        for ((c1, c2), sb) in sbsInPairComp.items2d():
            if geneNameInvolvedInAMistake == sb.la[0][0]:
                positionsOnGenomes.append((myLightGenomes.GeneP(c1, sb.l1[0][0]), myLightGenomes.GeneP(c2, sb.l2[0][0])))
            if geneNameInvolvedInAMistake == sb.la[-1][0]:
                positionsOnGenomes.append((myLightGenomes.GeneP(c1, sb.l1[-1][0]), myLightGenomes.GeneP(c2, sb.l2[-1][0])))
    elif mistake == 'Fn':
        posSbsGenome = simulatedSbsGenome.getPosition(gn, default=None)
        assert posSbsGenome is not None
        # assert that the gene is at the extremity of a sb of reference (chromExtremity mode)
        assert posSbsGenome.idx in {0, len(simulatedSbsGenome[posSbsGenome.c]) - 1}
        positionsOnGenomes.append((genome1.getPosition(gn, default=None), genome2.getPosition(gn, default=None)))
        # print >> sys.stderr, "%s on G1:%s and G2:%s" % (mistake, str(possGenome1), str(possGenome2))
    return positionsOnGenomes

    # distribLenSbsWithFn = collections.defaultdict(int)
    # for (geneN, _) in sFn:
    #     geneP = sbsGenome.getPosition(geneN, default=None)
    #     if geneP:
    #         lenSbFn = len(sbsGenome[geneP.c])
    #         print >> sys.stderr, "Fp: gene %s in chrR %s of len=%s at idx %s" % (geneN, geneP.c, lenSbFn, geneP.idx)
    #         distribLenSbsWithFn[lenSbFn] += 1
    # sumLens = sum(distribLenSbsWithFn.values())
    # percentageLen1 = 100 * float(distribLenSbsWithFn[1]) / sumLens if sumLens != 0 else 'None'
    # print >> sys.stderr, "(%s%% of len1), distribLenSbsWithFn=%s" %\
    #                      (percentageLen1,
    #                       sorted([(l, nb) for (l, nb) in distribLenSbsWithFn.iteritems()], key=lambda x: x[0]))
    # return distribLenSbsWithFn

def editSbs(referenceSbs, sbs, removeSbsOfLen1=False, reduceToSameGeneNames=False):
    if removeSbsOfLen1:
        sbs.removeChrsStrictlySmallerThan(2)
        assert all(len(sb) > 1 for sb in sbs.values())
        assert all(len(sb) > 1 for sb in referenceSbs.values())

    if reduceToSameGeneNames:
        sgSbs = sbs.getGeneNames(checkNoDuplicates=False)
        sgReferenceSbs = referenceSbs.getGeneNames()
        # reduce to the same set of genes
        setGenesInBoth = sgSbs & sgReferenceSbs
        setGenesToRemove = (sgSbs | sgReferenceSbs) - setGenesInBoth
        (sbs, _, (nbRemovedSbs1, nbRemovedGenes1)) =\
            myMapping.remapFilterGeneContent(sbs, setGenesToRemove)
        (nbRemovedSbs1, nbRemovedGenes1)
        (referenceSbs2, _, (nbRemovedSbs2, nbRemovedGenes2)) =\
             myMapping.remapFilterGeneContent(referenceSbs, setGenesToRemove)
        print >> sys.stderr, (nbRemovedSbs2, nbRemovedGenes2)
    else:
        referenceSbs2 = referenceSbs
    return (referenceSbs2, sbs)

def analyseFamiliesInSyntenyBlocks(sbsGenome, sbsGenomeTrue):
    assert isinstance(sbsGenomeTrue, myLightGenomes.LightGenome)
    assert isinstance(sbsGenome, myLightGenomes.LightGenome)
    setFamiliesTrue = sbsGenomeTrue.getGeneNames(asA=set, checkNoDuplicates=True)
    listFamilyNamesOfGenes = sbsGenome.getGeneNames(asA=list, checkNoDuplicates=False)
    setFamilies = set(listFamilyNamesOfGenes)
    nbOfDifferentFamilies = len(setFamilies)

    countFamilyNames = collections.Counter(listFamilyNamesOfGenes)
    assert nbOfDifferentFamilies == len(countFamilyNames.keys())
    assert setFamilies == set(countFamilyNames.keys())

    #sFp = setFamilies - setFamiliesTrue
    sFn = setFamiliesTrue - setFamilies
    # sTp = setFamiliesTrue & setFamilies
    # print >> sys.stderr, 'families in sbs but not in ref sbs = %s' % len(setFamilyNamesNotInRefSbs)
    # print >> sys.stderr, 'families in ref sbs but not in sbs = %s' % len(sFn)
    # print >> sys.stderr, 'families in ref sbs and in sbs = %s' % len(setCoreFamilies)

    # for ideal sbs, each family of the anc genes in ref sbs should be in one copy in sbs
    # nbSurplusCopiesOfFamilyName = sum(countFamiliesIdentified.values()) - nbOfDifferentFamiliesIdentified
    # if nbOfDifferentFamilies != 0 :
    #     percentageSurplusCopiesOfFamilyName = 100 * (float(nbSurplusCopiesOfFamilyName) / float(nbOfDifferentFamilies))
    # else:
    #     percentageSurplusCopiesOfFamilyName = 100
    nbSurplusCopiesOfFamilyInstances = sum([v for fn, v in countFamilyNames.iteritems() if fn in setFamiliesTrue]) - len(set(countFamilyNames.keys()) & setFamiliesTrue)
    print >> sys.stderr, 'nb surplus copies of (family names in ref sbs) = %s' % nbSurplusCopiesOfFamilyInstances

    # for ideal sbs, only gene families of ref sbs should be in sbs
    nbGenesWithNoCorrespondingTrueFamily = len([g_fn for g_fn in listFamilyNamesOfGenes if g_fn not in setFamiliesTrue])
    print >> sys.stderr, 'nb genes with no corresponding family name in ref sbs = %s' % nbGenesWithNoCorrespondingTrueFamily

    # genes in sbs + missing anc genes - (surplus copies of anc genes in ref sbs + anc genes not in ref sbs) = anc genes in ref sbs
    assert len(listFamilyNamesOfGenes) + len(sFn) - \
           (nbSurplusCopiesOfFamilyInstances + nbGenesWithNoCorrespondingTrueFamily) \
           == len(setFamiliesTrue)

# example
# variableArg = (('imr', ('imcs', 'gmmi'), 'tm'), [(False, (False, None), None),
#                                                 (True, (False, None), None),
#                                                 (False, (True, 0), None),
#                                                 (True, (True, 0), None),
#                                                 (True, (True, 1), None),
#                                                 (True, (True, 1), 4)])

# variableArg2 is a synthetic version of variableArg
# identify micro-rearrangements (imr)
# identify monogenic conserved segments (imcs) (either withing gaps of synteny blocks or at extremities of sbs)
# truncation
# variableArg2 = (('imr', 'imcs', 't'), [('-', '-', '-'),
#                                        ('+', '-', '-'),
#                                        ('-', '0', '-'),
#                                        ('+', '0', '-'),
#                                        ('+', '1', '-'),
#                                        ('+', '1', '4')])

def variableIntoStr(v):
    v_str = None
    if v is True:
        v_str = '+'
    elif v is False:
        v_str = '-'
    elif v is None:
        v_str = '-'
    elif isinstance(v, int):
        v_str = str(v)
    else:
        raise ValueError('v = %s is not a valid value' % v)
    return v_str

def variableArgIntoVariableArgString(variableArg):
    template = variableArg[0]
    vss = variableArg[1]
    assert all([len(vs) == len(template) for vs in vss])
    vss_str = []
    for vs in vss:
        (imr, (imcs, gmmi), tm) = vs
        new_vs = [None,None,None]
        new_vs[0] = imr
        new_vs[1] = gmmi if imcs is True else False
        new_vs[2] = tm
        vs_str = []
        for v in new_vs:
            vs_str.append(variableIntoStr(v))
        assert len(vs_str) == len(template)
        vss_str.append(vs_str)
    template = ('imr', 'imcs', 't')
    return (template, vss_str)


def ensureNamesDifferBetweenGenome1AndGenome2(genome1, genome2, families):
    # be sure that the names of genes differ between the two genomes
    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    newGenome1 = myLightGenomes.LightGenome(genome1)
    newGenome2 = myLightGenomes.LightGenome(genome2)
    for chr in newGenome1.keys():
        newGenome1[chr] = [myLightGenomes.OGene('1' + g.n, g.s) for g in genome1[chr]]
    for chr in newGenome2.keys():
        newGenome2[chr] = [myLightGenomes.OGene('2' + g.n, g.s) for g in genome2[chr]]
    newFamilies = myLightGenomes.Families()
    for family in families:
        newDns = set()
        for dn in family.dns:
            newDns.add('1' + dn)
            newDns.add('2' + dn)
        newFamilies.addFamily(myLightGenomes.Family(family.fn, newDns))
    return (newGenome1, newGenome2, newFamilies)