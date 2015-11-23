#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import sys

from utils import myPhylTree, myTools, myDiags, myLightGenomes, myIntervals


__doc__ = """ Extract ancGenes of sbs extremities """
__author__ = 'jlucas'

modesOrthos = list(myDiags.FilterType._keys)
arguments = myTools.checkArgs(
    # allPairwiseComparisons and for each comparison the ancGene file for gene
    # family definition.
    [('speciesTree', file)],
    [('genomes', str, 'genes.%s.list.bz2'),
     ('ancGenes', str, 'ancGenes.%s.list.bz2'),
     ("filterType", str, 'InBothGenomes'),
     ("tandemGapMax", int, 0),
     ("gapMax", str, 'None'),
     ("distinguishMonoGenicDiags", bool, True),
     ('distanceMetric', str, 'CD'),
     ('pThreshold', str, 'None'),
     ('gapMaxMicroInv', str, '0'),
     ('identifyMonoGenicInvs', bool, False),
     ('identifyMicroRearrangements', bool, True),
     ('truncationMax', str, 'None'),
     ("minChromLength", int, 2),
     ("sameStrand", bool, True),
     ('nbHpsRecommendedGap', int, 2), ('targetProbaRecommendedGap', float, 0.01),
     ('validateImpossToCalc_mThreshold', int, 3),
     ("out:sbsExtremities", str, 'res/sbsExtremities.%s.%s.%s.list'),
     ("out:sbs", str, 'res/sbs.%s.%s.list.bz2'),
     ('verbose', bool, True)
     ],
    __doc__
    )

for (argN, tpe) in [('gapMax', int), ('truncationMax', int), ('gapMaxMicroInv', int), ('pThreshold', float)]:
    if arguments[argN] == 'None':
        arguments[argN] = None
    else:
        try:
            arguments[argN] = tpe(arguments[argN])
        except:
            raise TypeError('%s is either an int or None' % argN)

filterType = myDiags.FilterType[modesOrthos.index(arguments["filterType"])]

# Define all the pairwise comparisons
allPComps = []
setOfExtantSpecies = set([])
phylTree = myPhylTree.PhylogeneticTree(arguments["speciesTree"])
for anc in phylTree.listAncestr:
    for (f1, f2) in itertools.combinations([f for (f, _) in phylTree.items[anc]], 2):
        l1 = phylTree.species[f1]
        l2 = phylTree.species[f2]
        for (sp1, sp2) in itertools.product(l1, l2):
            setOfExtantSpecies = setOfExtantSpecies | {sp1, sp2}
            allPComps.append((sp1, sp2, anc))
nbExtantSpecies = len(setOfExtantSpecies)

OGene = myLightGenomes.OGene
allSbsExtremities = {}
for (sp1, sp2, anc) in allPComps:

    genome1 = myLightGenomes.LightGenome(arguments['genomes'] % sp1)
    genome2 = myLightGenomes.LightGenome(arguments['genomes'] % sp2)
    ancGenes = myLightGenomes.Families(arguments['ancGenes'] % anc)

    allSbsExtremities[(sp1, sp2, anc)] = set()

    sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, ancGenes,
        filterType=filterType,
        tandemGapMax=arguments['tandemGapMax'],
        gapMax=arguments["gapMax"],
        distinguishMonoGenicDiags=arguments["distinguishMonoGenicDiags"],
        distanceMetric=arguments['distanceMetric'],
        pThreshold=arguments['pThreshold'],
        gapMaxMicroInv=arguments["gapMaxMicroInv"],
        identifyMonoGenicInvs=arguments["identifyMonoGenicInvs"],
        identifyMicroRearrangements=arguments['identifyMicroRearrangements'],
        truncationMax=arguments['truncationMax'],
        minChromLength=arguments["minChromLength"],
        sameStrand=arguments["sameStrand"],
        nbHpsRecommendedGap=arguments['nbHpsRecommendedGap'],
        targetProbaRecommendedGap=arguments['targetProbaRecommendedGap'],
        validateImpossToCalc_mThreshold=arguments['validateImpossToCalc_mThreshold'],
        verbose=arguments['verbose'])

    # record sbs
    myDiags.printSbsFile(sbsInPairComp, genome1, genome2, stream=open(arguments['out:sbs'] % (sp1, sp2), 'w'))

    for (_, sb) in sbsInPairComp.items2d():
        # oriented gene at the left-extremity of the sb
        ogl = OGene(sb.la[0][0], sb.la[0][1])
        # oriented gene at the right-extremity of the sb
        ogr = OGene(sb.la[-1][0], sb.la[-1][1])
        # gene extremity of the left-end of the sb
        gel = myIntervals.geneExtremityFromGene(ogl, -1)
        allSbsExtremities[(sp1, sp2, anc)].add(gel)
        # gene extremity of the right-end of the sb
        ger = myIntervals.geneExtremityFromGene(ogr, +1)
        allSbsExtremities[(sp1, sp2, anc)].add(ger)

for ((sp1, sp2, anc), sbsExtremities) in allSbsExtremities.iteritems():
    with open(arguments['out:sbsExtremities'] % (sp1, sp2, anc), 'w') as f:
        for (gn, hort) in sbsExtremities:
            # hort : head or tail
            # hort == t if tail
            # hort == d if head
            print >> f, "%s\t%s" % (gn, -1 if hort == 't' else +1)