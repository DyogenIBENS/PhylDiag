#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

__doc__ = """
Wrapper for PhylDiag Library
"""

import sys

import utils.myTools as myTools
import utils.myDiags as myDiags
import utils.myGenomes as myGenomes

# Arguments
modesOrthos = list(myDiags.FilterType._keys)
arguments = myTools.checkArgs( \
                [("genome1",file), ("genome2",file), ("ancGenes",file)], \
                [("tandemGapMax", int, 0),
                 ("gapMax",str,'None'),
                 ("sameStrand",bool,True), ("filterType",str,modesOrthos), ("minChromLength",int,1), ('distanceMetric',str,'CD'), ('pThreshold',float, 0.001),\
                ('nbHpsRecommendedGap',int,2), ('targetProbaRecommendedGap',float,0.01),\
                ('validateImpossToCalc_mThreshold',int,3),\
                ('multiProcess',bool,True),\
                ('verbose',bool,False)], \
                #, ("computeDiagsWithoutGenesOnlyImplyedInDiagsOfLengthSmallerOrEqualTo",int,-1)], \
                __doc__ \
                )

if arguments['gapMax'] == 'None':
    arguments['gapMax']=None
else:
    try:
        arguments['gapMax']=int(arguments['gapMax'])
    except:
        raise TypeError('gapMax is either an int or None')

genome1 = myGenomes.Genome(arguments["genome1"])
print >> sys.stderr, "Genome1"
print >> sys.stderr, "nb of Chr = ", len(genome1.chrList[myGenomes.ContigType.Chromosome])," nb of scaffolds = ", len(genome1.chrList[myGenomes.ContigType.Scaffold])
genome2 = myGenomes.Genome(arguments["genome2"])
print >> sys.stderr, "Genome2"
print >> sys.stderr, "nb of Chr = ", len(genome2.chrList[myGenomes.ContigType.Chromosome])," nb of scaffolds = ", len(genome2.chrList[myGenomes.ContigType.Scaffold])
ancGenes = myGenomes.Genome(arguments["ancGenes"])
filterType = myDiags.FilterType[modesOrthos.index(arguments["filterType"])]
statsDiags = []

print >> sys.stderr, "Begining of the extraction of synteny blocks"
listOfSbs =\
    list(myDiags.extractSbsInPairCompGenomes(genome1, genome2, ancGenes,
                                             tandemGapMax=arguments['tandemGapMax'],
                                             gapMax=arguments["gapMax"],
                                             distanceMetric=arguments['distanceMetric'],
                                             pThreshold=arguments['pThreshold'],
                                             filterType=filterType,
                                             minChromLength=arguments["minChromLength"],
                                             consistentSwDType=arguments["sameStrand"],
                                             nbHpsRecommendedGap=arguments['nbHpsRecommendedGap'],
                                             targetProbaRecommendedGap=arguments['targetProbaRecommendedGap'],
                                             validateImpossToCalc_mThreshold=arguments['validateImpossToCalc_mThreshold'],
                                             multiProcess=arguments['multiProcess'],
                                             verbose=arguments['verbose']))
print >> sys.stderr, "End of the synteny block research"
listOfSbs.sort(key=lambda x: len(x[2]), reverse=True)

myDiags.printSbsFile(listOfSbs, genome1, genome2)
