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

import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myDiags

# Arguments
modesOrthos = list(utils.myDiags.FilterType._keys)
arguments = utils.myTools.checkArgs( \
		[("genome1",file), ("genome2",file), ("ancGenes",file)], \
		[("gapMax",str,'None'), ("sameStrand",bool,True), ("filterType",str,modesOrthos), ("minChromLength",int,1), ('distanceMetric',str,'CD'), ('pThreshold',float, 0.001),\
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

genome1 = utils.myGenomes.Genome(arguments["genome1"])
print >> sys.stderr, "Genome1"
print >> sys.stderr, "nb of Chr = ", len(genome1.chrList[utils.myGenomes.ContigType.Chromosome])," nb of scaffolds = ", len(genome1.chrList[utils.myGenomes.ContigType.Scaffold])
genome2 = utils.myGenomes.Genome(arguments["genome2"])
print >> sys.stderr, "Genome2"
print >> sys.stderr, "nb of Chr = ", len(genome2.chrList[utils.myGenomes.ContigType.Chromosome])," nb of scaffolds = ", len(genome2.chrList[utils.myGenomes.ContigType.Scaffold])
ancGenes = utils.myGenomes.Genome(arguments["ancGenes"])
filterType = utils.myDiags.FilterType[modesOrthos.index(arguments["filterType"])]
statsDiags = []

print >> sys.stderr, "Begining of the extraction of synteny blocks"
listOfDiags = list(utils.myDiags.extractSbsInPairCompGenomes(genome1, genome2, ancGenes, gapMax=arguments["gapMax"], distanceMetric=arguments['distanceMetric'], pThreshold=arguments['pThreshold'], filterType=filterType, minChromLength=arguments["minChromLength"], consistentSwDType=arguments["sameStrand"], nbHpsRecommendedGap=arguments['nbHpsRecommendedGap'], targetProbaRecommendedGap=arguments['targetProbaRecommendedGap'], validateImpossToCalc_mThreshold=arguments['validateImpossToCalc_mThreshold'], multiProcess=arguments['multiProcess'], verbose=arguments['verbose']))
print >> sys.stderr, "End of the synteny block research"
listOfDiags.sort(key=lambda x: len(x[2]))
lenListOfDiags = len(listOfDiags)

for ((c1,d1),(c2,d2),daa,pVal) in listOfDiags:

	l = len(daa)
	if l < arguments["minChromLength"]:
		continue
	statsDiags.append(l)

	res = [l, \
		c1," ".join(genome1.lstGenes[c1][i1].names[0] for (i1,_) in d1), \
		c2," ".join(genome2.lstGenes[c2][i2].names[0] for (i2,_) in d2) ]
	# The ancestral orientation, width and high of the homology packs are removed to fit Matthieu's format
	res.append(utils.myFile.myTSV.printLine([a[0] for a in daa], " "))
	res.append(utils.myFile.myTSV.printLine([a[2] for a in daa], " "))
	res.append(utils.myFile.myTSV.printLine([a[3] for a in daa], " "))

	if arguments["sameStrand"]:
		# orientation of the hp sign
		res.append(utils.myFile.myTSV.printLine([a[1] for a in daa], " "))

	# We add the pValue
	res.append(utils.myFile.myTSV.printLine([pVal], " "))

	print utils.myFile.myTSV.printLine(res)

print >> sys.stderr, "Distribution of the synteny block lengths", utils.myMaths.myStats.txtSummary(statsDiags)
