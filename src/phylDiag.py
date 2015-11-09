#!/usr/bin/python
# -*- coding: utf-8 -*-

# PhylDiag version 2.0 (6/11/2015)
# python 2.7
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

import collections

__doc__ = """
Wrapper for the PhylDiag Library in LibsDyogen
"""

import sys

import utils.myTools as myTools
import utils.myDiags as myDiags
import utils.myLightGenomes as myLightGenomes

# Arguments
arguments = myTools.checkArgs(
    [("genome1", file),
     ("genome2", file),
     ("families", file)],
    myDiags.defaultArgsPhylDiag + [('removeUnofficialChromosomes', bool, True)],
    #, ("computeDiagsWithoutGenesOnlyImplyedInDiagsOfLengthSmallerOrEqualTo",int,-1)], \
    __doc__
)
kwargs = myDiags.defaultKwargsPhylDiag(arguments=arguments)
kwargs['verbose'] = True

genome1 = myLightGenomes.LightGenome(arguments["genome1"])
if arguments['removeUnofficialChromosomes']:
    genome1.removeUnofficialChromosomes()
print >> sys.stderr, "Genome1"
print >> sys.stderr, "Nb of Chr = ", len(genome1.keys())
genome2 = myLightGenomes.LightGenome(arguments["genome2"])
if arguments['removeUnofficialChromosomes']:
    genome2.removeUnofficialChromosomes()
print >> sys.stderr, "Genome2"
print >> sys.stderr, "Nb of Chr = ", len(genome2.keys())
families = myLightGenomes.Families(arguments["families"])
#filterType = myDiags.FilterType[list(myDiags.FilterType._keys).index(arguments["filterType"])]
statsDiags = []

print >> sys.stderr, "Beginning of the extraction of synteny blocks"
sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, families, **kwargs)
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