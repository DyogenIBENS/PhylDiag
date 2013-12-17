#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.0
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work
# under the terms of the GNU General Public License, version 2 or later

__doc__ = """
	Lit les arbres de proteines et cree les fichiers de genes ancestraux
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
#import utils.myProteinTree

import imp
utils.myProteinTree = imp.load_source('utils.myProteinTree', '/workspace4/jlucas/Workflow/Scripts/Updated/myProteinTree.py')


# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("proteinTree",file)], [("OUT.ancGenesFiles",str,""), ("reuseNames",bool,False)], __doc__ )

#phylTree = speciesTree
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

#duplication counter (using a collection and a default dict of ints allow us to avoid creating a key and a first element. Using .expand, even for the first time the key is used will create an entry key. 
dupCount = collections.defaultdict(int)
def futureName(name, dup):
	if dup >= 2:
		dupCount[name] += 1
		#If there is a duplication we need to add a suffix
		return name + utils.myProteinTree.getDupSuffix(dupCount[name], False)
	else:
		return name

#################################################
# Retrouve les vraies racines dans les familles #
#################################################
def getRoots(node, previousAnc, lastWrittenAnc):  

	newAnc = tree.info[node]['taxon_name'] 
	
	# Indique les noms des ancetres depuis le dernier lu, ici on ne garde que le newLastWritten
	(_,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, tree.info[node]['Duplication'] >=2)

	if isroot:
		return [node]

	# Les genes des descendants
	subRoots = []
	for (g,_) in tree.data.get(node,[]):
		subRoots.extend( getRoots(g, newAnc, newLastWritten) )
	return subRoots

count = collections.defaultdict(int)

###########################################
# Sauvegarde toutes les familles de genes #
###########################################
def extractGeneFamilies(node, baseName, previousAnc, lastWrittenAnc):

	newAnc = tree.info[node]['taxon_name']

	(toWrite,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, tree.info[node]['Duplication'] >= 2)

	if isroot and (previousAnc != None): # previousAnc != None evite la racine de l'arbre de genes. isroot est vrai pour les noeuds de speciation du premier ancetre, ou les noeuds de duplication juste apres le 1er ancetre  
		if not arguments["reuseNames"]:
			baseName = baseName.split(".")[0]
		count[baseName] += 1
		currName = baseName + utils.myProteinTree.getDupSuffix(count[baseName], True)
	else:
		currName = baseName
	tree.info[node]['family_name'] = currName # ! on modifie 'family name' de tree alors qu'il n'est meme pas en parametre

	# Les genes des descendants
	if node in tree.data: # vrai si le noeud est dans l'arbre des proteines (si le noeud est une feuille il n'est pas dans tree.data)
		allGenes = []
		for (g,_) in tree.data[node]: # {(g.a,len_a),(g_b,len_b),...} 
			allGenes.extend( extractGeneFamilies(g, futureName(currName, tree.info[node]['Duplication']), newAnc, newLastWritten) )

	else: # Quand c'est une feuille!
		allGenes = [ tree.info[node]["gene_name"] ]

	for a in toWrite: # 'a'= Nom de l'ancetre a ecrire 
	 	geneFamilies[a].append( [currName] + allGenes ) # Ici on ecrit le nom du gene de l'Anc suivi des noms des genes des especes descendantes (et ca pour toutes les especes dans toWrite)
		#! geneFamilies est definie dans le main, on le modifie alors qu'il n'est pas en parametre 

	return allGenes # utile pour la definition par recurrence

geneFamilies = collections.defaultdict(list)
for tree in utils.myProteinTree.loadTree(arguments["proteinTree"]): # pour tous les arbres dans la foret d'arbre de genes
	extractGeneFamilies(tree.root, tree.info[tree.root]["tree_name"], None, None) # cette fonction modifie tree et geneFamilies bien que tree et geneFamilies ne sont pas des parametre
	tree.printTree(sys.stdout)

for (anc,lst) in geneFamilies.iteritems():
	print >> sys.stderr, "Ecriture des familles de %s ..." % anc,
	f = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
	for gg in lst:
		print >> f, " ".join(gg)
	f.close()
	print >> sys.stderr, len(lst), "OK"


