#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# This script was originally written by Matthieu Muffato, it has been integrated
# to the PhylDiag package for practical uses

__doc__ = """
        Plot the homology matrix between two genomes
"""

from utils import myTools, myLightGenomes, myGenomesDrawer, myDiags

arguments = myTools.checkArgs(
        [("genome1", file), ("genome2", file), ("families", file)],
        [('removeUnofficialChromosomes', bool, True),
         ('withSbs', bool, True),
         ('filterType', str, 'InBothGenomes')],
        __doc__)

filterType = list(myDiags.FilterType._keys)
filterType = myDiags.FilterType[filterType.index(arguments["filterType"])]

genome1 = myLightGenomes.LightGenome(arguments['genome1'])
genome2 = myLightGenomes.LightGenome(arguments['genome2'])
families = myLightGenomes.Families(arguments['families'])

if arguments['removeUnofficialChromosomes']:
    genome1.removeUnofficialChromosomes()
    genome2.removeUnofficialChromosomes()

#FIXME :
# write chromosomes names

sbsInPairComp = None
if arguments['withSbs']:
    sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, families,
                                                        filterType=filterType,
                                                        tandemGapMax=0,
                                                        gapMax=5,
                                                        distanceMetric='CD',
                                                        distinguishMonoGenicDiags=False,
                                                        pThreshold=1.0,
                                                        gapMaxMicroInv=1,
                                                        identifyMonoGenicInversion=False,
                                                        identifyBreakpointsWithinGaps=False,
                                                        overlapMax=20,
                                                        consistentSwDType=True,
                                                        minChromLength=2,
                                                        nbHpsRecommendedGap=2,
                                                        targetProbaRecommendedGap=0.01,
                                                        validateImpossToCalc_mThreshold=3,
                                                        optimisation='cython',
                                                        verbose=False)

myGenomesDrawer.wholeGenomeHomologyMatrices(genome1, genome2, families,
                                            inSbsInPairComp=sbsInPairComp,
                                            filterType=filterType,
                                            minChromLength=0,
                                            tandemGapMax=4,
                                            scaleFactorRectangles=10,
                                            outputFileName='res/pairwiseCompGenomesHMs.svg',
                                            maxWidth=100,
                                            maxHeight=100)