#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr or hrc@ens.fr
# Licences GLP v3 and CeCILL v2

import collections

__doc__ = """
Wrapper for the PhylDiag Library in LibsDyogen that uses myTools.checkArgs
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
    myDiags.defaultArgsPhylDiag,
    #, ("computeDiagsWithoutGenesOnlyImplyedInDiagsOfLengthSmallerOrEqualTo",int,-1)], \
    __doc__
)
kwargs = myDiags.defaultKwargsPhylDiag(arguments=arguments)
kwargs['verbose'] = True

print >> sys.stderr, 'Warning: remember that phylDiag do not removeUnofficialChromosomes by default'

genome1 = myLightGenomes.LightGenome(arguments["genome1"])
print >> sys.stderr, "Genome1"
print >> sys.stderr, "Nb of Chr = ", len(genome1.keys())
genome2 = myLightGenomes.LightGenome(arguments["genome2"])
print >> sys.stderr, "Genome2"
print >> sys.stderr, "Nb of Chr = ", len(genome2.keys())
families = myLightGenomes.Families(arguments["families"])

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