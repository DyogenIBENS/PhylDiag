#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France
import collections

__doc__ = """
Wrapper for PhylDiag Library
"""

import sys

import utils.myTools as myTools
import utils.myDiags as myDiags
#import utils.myGenomes as myGenomes
import utils.myLightGenomes as myLightGenomes

# Arguments
arguments = myTools.checkArgs(
                [("genome1", file),
                 ("genome2", file),
                 ("families", file)],
                [("filterType", str, 'InBothGenomes'),
                 ("tandemGapMax", int, 0),
                 ("gapMax", str, 'None'),
                 ("distinguishMonoGenicDiags", bool, True),
                 ('distanceMetric', str, 'CD'),
                 ('pThreshold', str, 'None'),
                 ('gapMaxMicroInv', str, '0'),
                 ('identifyBreakpointsWithinGaps', bool, True),
                 ('overlapMax', str, 'None'),
                 ("minChromLength", int, 2),
                 ("sameStrand", bool, True),
                 ('nbHpsRecommendedGap', int, 2), ('targetProbaRecommendedGap', float, 0.01),
                 ('validateImpossToCalc_mThreshold', int, 3),
                 # The multiprocess does not seem to work well for most of the data. It work well only for some data
                 # with ~ 800 contigs
                ('multiProcess', bool, False),
                ('verbose', bool, False)],
                #, ("computeDiagsWithoutGenesOnlyImplyedInDiagsOfLengthSmallerOrEqualTo",int,-1)], \
                __doc__
                )

for (argN, tpe) in [('gapMax', int), ('overlapMax', int), ('gapMaxMicroInv', int), ('pThreshold', float)]:
    if arguments[argN] == 'None':
        arguments[argN] = None
    else:
        try:
            arguments[argN] = tpe(arguments[argN])
        except:
            raise TypeError('%s is either an int or None' % argN)

genome1 = myLightGenomes.LightGenome(arguments["genome1"])
print >> sys.stderr, "Genome1"
print >> sys.stderr, "Nb of Chr = ", len(genome1.keys())
genome2 = myLightGenomes.LightGenome(arguments["genome2"])
print >> sys.stderr, "Genome2"
print >> sys.stderr, "Nb of Chr = ", len(genome2.keys())
families = myLightGenomes.Families(arguments["families"])
filterType = list(myDiags.FilterType._keys)
filterType = myDiags.FilterType[filterType.index(arguments["filterType"])]
statsDiags = []

print >> sys.stderr, "Begining of the extraction of synteny blocks"
sbsInPairComp = \
    myDiags.extractSbsInPairCompGenomes(genome1, genome2, families,
                                        filterType=filterType,
                                        tandemGapMax=arguments['tandemGapMax'],
                                        gapMax=arguments["gapMax"],
                                        distinguishMonoGenicDiags=arguments["distinguishMonoGenicDiags"],
                                        distanceMetric=arguments['distanceMetric'],
                                        pThreshold=arguments['pThreshold'],
                                        gapMaxMicroInv=arguments["gapMaxMicroInv"],
                                        identifyBreakpointsWithinGaps=arguments['identifyBreakpointsWithinGaps'],
                                        overlapMax=arguments['overlapMax'],
                                        minChromLength=arguments["minChromLength"],
                                        consistentSwDType=arguments["sameStrand"],
                                        nbHpsRecommendedGap=arguments['nbHpsRecommendedGap'],
                                        targetProbaRecommendedGap=arguments['targetProbaRecommendedGap'],
                                        validateImpossToCalc_mThreshold=arguments['validateImpossToCalc_mThreshold'],
                                        multiProcess=arguments['multiProcess'],
                                        verbose=arguments['verbose'])
print >> sys.stderr, "End of the synteny block research"
nbOfSbs = sum([len(sbsInPairComp[c1][c2]) for (c1, c2) in sbsInPairComp.keys2d()])
print >> sys.stderr, "Number of synteny blocks = %s" % nbOfSbs
if nbOfSbs > 0:
    sbLens = sorted([len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()])
    nbSbsOfLen = collections.Counter(sbLens)
    distribSbLens = [" %s:%s" % (nbSbsOfLen[sbLen], sbLen) for sbLen in sorted(nbSbsOfLen.keys())]
    distribSbLens = distribSbLens[:5] + ["..."] + distribSbLens[-3:]
    print >> sys.stderr, "Distribution of sb lengths (nbSbs:length) = %s" % ",".join(distribSbLens)

# FIXME: create an empty file or do not create any file when no sbs
myDiags.printSbsFile(sbsInPairComp, genome1, genome2, sortByDecrLengths=True)
