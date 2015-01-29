#!/usr/bin/python
# -*- coding: utf-8 -*-
# Wrapper for 'homologyteams'
# http://euler.slu.edu/~goldwasser/homologyteams/
# Identifying Conserved Gene Clusters in the Presence of Homology Families
# Xin He and Michael H. Goldwasser
# Journal of Computational Biology, 12(6-7):638-656, 2005.)

# * Project: homologyteams
# * Authors: Michael H. Goldwasser (goldwamh@slu.edu) and Xin He (xinhe2@uiuc.edu)
# * Version: 1.0 (May 2004)
# *
# * Copyright (C) 2004  Michael H. Goldwasser and Xin He
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Copyright of the wrapper
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

__doc__ = """
Wrapper for gene teams library 'homologyteams'
http://euler.slu.edu/~goldwasser/homologyteams/
Identifying Conserved Gene Clusters in the Presence of Homology Families
Xin He and Michael H. Goldwasser
Journal of Computational Biology, 12(6-7):638-656, 2005.)
"""

import sys

import utils.myTools as myTools
import utils.myGeneTeams as myGeneTeams
import utils.myLightGenomes as myLightGenomes

# Arguments
modesOrthos = list(myGeneTeams.FilterType._keys)
arguments = myTools.checkArgs(
    [("genome1", file),
     ("genome2", file),
     ("families", file)],
    [("tandemGapMax", int, 0),
     ("gapMax", str, 'None'),
     ("filterType", str, modesOrthos),
     ("minChromLength", int, 1),
     ('verbose', bool, False)],
    __doc__
    )

if arguments['gapMax'] == 'None':
    arguments['gapMax'] = None
else:
    try:
        arguments['gapMax'] = int(arguments['gapMax'])
    except:
        raise TypeError('gapMax is either an int or None')

genome1 = myLightGenomes.LightGenome(arguments["genome1"])
print >> sys.stderr, "Genome1"
print >> sys.stderr, "Nb of Chr = ", len(genome1)
genome2 = myLightGenomes.LightGenome(arguments["genome2"])
print >> sys.stderr, "Genome2"
print >> sys.stderr, "Nb of Chr = ", len(genome2)
families = myLightGenomes.Families(arguments["families"])
filterType = myGeneTeams.FilterType[modesOrthos.index(arguments["filterType"])]
statsDiags = []

print >> sys.stderr, "Begining of the extraction of gene teams"

# in order to test easily (faster execution) on the comparison of both X chromosomes
test = False
if test is True:
    genome1 = genome1.intoDict()
    genome2 = genome2.intoDict()
    g1 = {}
    g1['X'] = genome1['X']
    g2 = {}
    g2['X'] = genome2['X']
    genome1 = g1
    genome2 = g2

listOfGts =\
    list(myGeneTeams.extractGtsInPairCompGenomes(genome1, genome2, families,
                                                 tandemGapMax=arguments['tandemGapMax'],
                                                 gapMax=arguments["gapMax"],
                                                 filterType=filterType,
                                                 minChromLength=arguments["minChromLength"],
                                                 verbose=arguments['verbose']))
print >> sys.stderr, "End of the gene team research"
# sort the list of gene teams by decreasing length
listOfGts.sort(key=lambda x: len(x[0]), reverse=True)
myGeneTeams.printGtsFile(listOfGts, genome1, genome2, families)
