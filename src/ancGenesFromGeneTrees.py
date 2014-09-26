#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

__doc__ = """
        extract the ancGenes from the forest of gene trees
"""

import sys
import collections

import libs.utils.myFile as myFile
import libs.utils.myTools as myTools
import libs.utils.myPhylTree as myPhylTree
import libs.utils.myProteinTree as myProteinTree

# arguments
arguments = myTools.checkArgs( [("speciesTree",file), ("geneTreeForest",file)], [("out:ancGenes",str,""), ("reuseNames",bool,False)], __doc__ )

speciesTree = myPhylTree.PhylogeneticTree(arguments["speciesTree"])

# duplication counter
dupCount = collections.defaultdict(int)
def futureName(name, dup):
    if dup >= 2:
        dupCount[name] += 1
        # if there is a duplication we need to add a suffix
        return name + myProteinTree.getDupSuffix(dupCount[name], False)
    else:
        return name

# finds out the roots in gene families
def getRoots(node, previousAnc, lastWrittenAnc):

    newAnc = tree.info[node]['taxon_name']

    # names of the ancestors since the last read, here only newLastWritten is kept
    (_,newLastWritten,isroot) = myProteinTree.getIntermediateAnc(speciesTree, previousAnc, lastWrittenAnc, newAnc, tree.info[node]['Duplication'] >=2)

    if isroot:
        return [node]

    # genes of children
    subRoots = []
    for (g,_) in tree.data.get(node,[]):
        subRoots.extend( getRoots(g, newAnc, newLastWritten) )
    return subRoots

count = collections.defaultdict(int)

# record all gene families
def extractGeneFamilies(node, baseName, previousAnc, lastWrittenAnc):

    newAnc = tree.info[node]['taxon_name']

    (toWrite,newLastWritten,isroot) = myProteinTree.getIntermediateAnc(speciesTree, previousAnc, lastWrittenAnc, newAnc, tree.info[node]['Duplication'] >= 2)

    if isroot and (previousAnc != None): # previousAnc != None avoid the root of the gene tree. isroot is true for the speciation nodes of the first ancestor, or the duplication nodes just after the first ancestor
        if not arguments["reuseNames"]:
            baseName = baseName.split(".")[0]
        count[baseName] += 1
        currName = baseName + myProteinTree.getDupSuffix(count[baseName], True)
    else:
        currName = baseName
    tree.info[node]['family_name'] = currName #FIXME ! 'family name' of 'tree' is modified whereas it is not even one of the parameter of the function...

    # genes of children
    if node in tree.data: # true if the node is in the forest of gene trees (if the node is a leaf it is not in tree.data)
        allGenes = []
        for (g,_) in tree.data[node]: # {(g.a,len_a),(g_b,len_b),...}
            allGenes.extend( extractGeneFamilies(g, futureName(currName, tree.info[node]['Duplication']), newAnc, newLastWritten) )

    else: # when the node is a leaf
        allGenes = [ tree.info[node]["gene_name"] ]

    for a in toWrite: # 'a'= name of the ancestor to print
        geneFamilies[a].append( [currName] + allGenes ) # write the name of the gene of Anc followed by the names of the genes of the children (this is done for all the species in toWrite)
        #FIXME geneFamilies is defined in the main, it is modified whereas it is not even a parameter

    return allGenes # for the recurrence

geneFamilies = collections.defaultdict(list)
for tree in myProteinTree.loadTree(arguments["geneTreeForest"]): # for all gene trees in the forest
    extractGeneFamilies(tree.root, tree.info[tree.root]["tree_name"], None, None) # FIXME this function modifies tree and geneFamilies even if tree and gene families are not parameters
    tree.printTree(sys.stdout)

for (anc,lst) in geneFamilies.iteritems():
    print >> sys.stderr, "Write %s family ..." % anc,
    f = myFile.openFile(arguments["out:ancGenes"] % speciesTree.fileName[anc], "w")
    for gg in lst:
        print >> f, " ".join(gg)
    f.close()
    print >> sys.stderr, len(lst), "OK"
