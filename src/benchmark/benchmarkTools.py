#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import numpy
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

import copy
import itertools
from utils import myFile
from utils import myTools
# TODO this function is the same as in MagSimus/src/libs/myBreakpointsAnalyser.py
# Once MagSimus is published, only conserve the function in MagSimus/src/libs/myBreakpointsAnalyser.py and remove this
# function and the subfunctions
# and edit links/imports properly
# instead of import 'from libs import myBreakpointsAnalyser.py' that already contains this function
def computePairwiseSyntenyBlocksOfMagSimus(speciesTree, pSimGenomes, pAncGenes, pReferenceSbs):

    # same function as MagSimus/src/libs/myEvents.performFission
    # fission on chromosome 'c', inner intergene of index 'x'
    # x is in [1, len(c)-1]
    # c -> c=c[:x] & newC=c[x:]
    def performFission(genome, (c, x), updateGenomeDict=False, keepOriginalGenome=False):
        assert isinstance(genome, myLightGenomes.LightGenome)
        if x == 0 or x == len(genome[c]):
            # do not create a new chromosome, do nothing
            newGenome = genome
        else:
            if keepOriginalGenome:
                newGenome = copy.copy(genome)
                newGenome[c] = copy.copy(genome[c])
            else:
                newGenome = genome
            newC = myLightGenomes.newChromName(genome)
            newGenome[newC] = newGenome[c][x:]
            del newGenome[c][x:]
            if updateGenomeDict:
                assert newGenome.withDict
                for (idx, g) in enumerate(newGenome[newC]):
                    newGenome.g2p[g.n] = myLightGenomes.GeneP(newC, idx)
        return newGenome

    # same function as MagSimus/src/libs/myEvents.combineTwoSbs
    def combineTwoSbs(sbs_1, sbs_2, verbose=False):
        """
        Combine two empirical sbs into one synthetic empirical sb

        :param sbs_1: sbs with positional orthologous gene names
        :param sbs_2: sbs with positional orthologous gene names
        :param verbose:
        :return: synthetic sbs
        """
        assert isinstance(sbs_1, myLightGenomes.LightGenome) and isinstance(sbs_2, myLightGenomes.LightGenome)
        sbs_1Gns = sbs_1.getGeneNames()
        sbs_2Gns = sbs_2.getGeneNames()

        ######################################
        # 1) reduce to the same set of genes #
        ######################################
        setGenesInBothSpecies = sbs_1Gns & sbs_2Gns
        setGenesToRemove = (sbs_1Gns | sbs_2Gns) - setGenesInBothSpecies
        (sbs_1, _, (nbRemovedSbsAnc, nbRemovedGenesAnc)) = myMapping.remapFilterGeneContent(sbs_1, setGenesToRemove)
        (sbs_2, _, (nbRemovedSbsDes, nbRemovedGenesDes)) = myMapping.remapFilterGeneContent(sbs_2, setGenesToRemove)
        assert sbs_1.getGeneNames() == sbs_2.getGeneNames()
        nbSbsRmed = nbRemovedSbsAnc + nbRemovedSbsDes

        ########################
        # 2) share breakpoints #
        ########################
        sbsExtremities1 = myIntervals.analyseGenomeIntoChromExtremities(sbs_1)
        sbsExtremities2 = myIntervals.analyseGenomeIntoChromExtremities(sbs_2)
        sbsExtremities = sbsExtremities1 | sbsExtremities2
        syntheticSbs = copy.deepcopy(sbs_1)
        assert all(len(sb) > 0 for sb in syntheticSbs.values())
        syntheticSbs.computeDictG2P()
        for breakpointGeneExtr in sbsExtremities:
            (c, x) = myIntervals.intergeneFromGeneExtremity(breakpointGeneExtr, syntheticSbs)
            syntheticSbs = performFission(syntheticSbs, (c, x), updateGenomeDict=True)
        assert all(len(sb) > 0 for sb in syntheticSbs.values())
        assert isinstance(syntheticSbs, myLightGenomes.LightGenome)

        return (syntheticSbs, nbSbsRmed)

    # same function as MagSimus/src/libs/myEvents.integrateAllBreakpoints
    def integrateAllBreakpoints(evolutivePathFromLCAtoD, empiricalSbs):
        """
        :param evolutivePathFromLCAtoD: example ['Amniota', 'Theria', 'Boreoeutheria', 'Euarchontoglires', 'Homo sapiens']
        :param empiricalSbs: dict of sbs with keys (Amniota, Theria, ..., 'Homo sapiens', ...)
        :param ancGenes: dict of ancGenes Families with keys (Amniota, Theria, ..., ...)
        :return: sbs
        """
        # ancestor -> descendant chain
        a_d_chain = myTools.myIterator.slidingTuple(evolutivePathFromLCAtoD, width=2)
        (lca, d) = a_d_chain.next()
        sbs = empiricalSbs[(lca, d)]
        assert all(len(sb) > 0 for sb in sbs.values())
        old_d = d
        # a: ancestor
        # d: descendant
        nbSbsRmed = 0
        for (a, d) in a_d_chain:
            assert a == old_d
            new_sbs = empiricalSbs[(a, d)]
            assert all(len(sb) > 0 for sb in new_sbs.values())
            # print >> sys.stdout, "in integrateAllBreakpoints %s->%s" % (a, d)
            (sbs, tmp_nbSbsRmed) = combineTwoSbs(sbs, new_sbs, verbose=True)
            nbSbsRmed += tmp_nbSbsRmed
            assert all(len(sb) > 0 for sb in sbs.values())
            old_d = d
        return (sbs, nbSbsRmed)

    # load families
    ancGenesOf = {}
    for s in speciesTree.listAncestr:
        ancGenesOf[s] = myLightGenomes.Families(pAncGenes % s)
    # load genomes
    genomeOf = {}
    for s in speciesTree.allNames:
        print >> sys.stderr, pSimGenomes % str(s)
        genomeOf[s] = myLightGenomes.LightGenome(pSimGenomes % str(s))
    # load sbs
    simSbsOf = {}

    def loadSbs(speciesTree, pSbs, parent):
        for (child, bLength) in speciesTree.items[parent]:
            simSbsOf[(parent, child)] = myLightGenomes.LightGenome(pSbs % (parent, str(child)))
            if child in speciesTree.items:
                # if child is an internal node of the species tree
                loadSbs(speciesTree, pSbs, child)
    loadSbs(speciesTree, pReferenceSbs, speciesTree.root)

    # for each pairwise comparison of two extant species compute the corresponding sbs
    for (sp1, sp2) in itertools.combinations(speciesTree.listSpecies, 2):
        lca = speciesTree.lastCommonAncestor([sp1, sp2])
        # print lca
        lca_genome = genomeOf[lca]
        nbGs_ini = len(lca_genome.getGeneNames(checkNoDuplicates=True))

        # sbs in the lineage from lca to sp1
        (sbs1, nbSbsRmed1) = integrateAllBreakpoints(speciesTree.dicLinks[lca][sp1], simSbsOf)
        # sbs in the lineage from lca to sp2
        (sbs2, nbSbsRmed2) = integrateAllBreakpoints(speciesTree.dicLinks[lca][sp2], simSbsOf)
        assert all(len(chrom) > 0 for chrom in sbs1.values())
        assert all(len(chrom) > 0 for chrom in sbs2.values())
        assert len(sbs1.getGeneNames(asA=list)) == len(sbs1.getGeneNames(asA=set))
        assert len(sbs2.getGeneNames(asA=list)) == len(sbs2.getGeneNames(asA=set))
        (sbs, nbSbsRmed3)  = combineTwoSbs(sbs1, sbs2, verbose=False)
        nbGsRmed = nbGs_ini - len(sbs.getGeneNames(checkNoDuplicates=True))
        nbSbsRmed = nbSbsRmed1 + nbSbsRmed2 + nbSbsRmed3
        # to have the species combination in alphabetical order
        speciesNames = sorted([sp1, sp2])
        print >> sys.stderr, "Computation of %s: %s sbs removed, %s genes removed" % \
                             (pReferenceSbs % (speciesNames[0], speciesNames[1]), nbSbsRmed, nbGsRmed)
        # sort sbs by decreasing sizes
        sbs.sort()
        with myFile.openFile(pReferenceSbs % (speciesNames[0], speciesNames[1]), 'w') as f:
            sbs.printIn(f)

Mylimits = numpy.zeros((3,3), dtype=object)
# very tight limits
# Mylimits[0][0] = (0.1,0.8)
# Mylimits[1][0] = (0.45,1)
# Mylimits[2][0] = (0.25,1)
# Mylimits[0][1] = (0.5,1)
# Mylimits[1][1] = (0.79,0.92)
# Mylimits[2][1] = (0.60,0.95)
# Mylimits[0][2] = (0.5,1)
# Mylimits[1][2] = (0.955,0.98)
# Mylimits[2][2] = (0.70,1)

Mylimits[0][0] = (0.1,1)
Mylimits[1][0] = (0.4,1)
Mylimits[2][0] = (0.2,1)

Mylimits[0][1] = (0.5,1)
Mylimits[1][1] = (0.7,1)
Mylimits[2][1] = (0.6,1)

Mylimits[0][2] = (0.5,1)
Mylimits[1][2] = (0.9,1)
Mylimits[2][2] = (0.7,1)