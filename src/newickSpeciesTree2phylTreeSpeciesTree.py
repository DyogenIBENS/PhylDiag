#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

""" Convert species tree (newick to phylTree) or (phylTree to newick) """

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs([("phylTree.conf",file)], [("fromNewick",bool,True)], "Convertit un arbre en phylogenetique depuis/vers le format Newick")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


if arguments["fromNewick"]:

	# Impression sous mon format, avec des indentations
	def do(node, indent):
		node = node.replace("*", "")
		node = node.replace("."," ")
		names = utils.myFile.myTSV.printLine([node] + [x for x in phylTree.commonNames.get(node,"") if isinstance(x, str) and (x != node)], delim="|")
		print ("\t" * indent) + "%s" % names
		if node in phylTree.items:
			for (f,_) in phylTree.items[node]:
				do(f, indent+1)

	do(phylTree.root, 0)

else:
	# Renvoie l'arbre au format Newick
	def convertToFlatFile(anc):

		a = phylTree.fileName[anc] # anc.replace(' ', '.')
		if anc in phylTree.listSpecies:
			return a
		else:
			return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for (e,l) in phylTree.items[anc]]) + ")%s|%d" % (a,phylTree.ages[anc])
	print convertToFlatFile(phylTree.root), ";"

