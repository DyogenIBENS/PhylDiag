# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import sys
import itertools
import collections

import myFile
import myTools

GeneSpeciesPosition = collections.namedtuple("GeneSpeciesPosition", ['species', 'chromosome', 'index'])


# class mother of all types of trees
class PhylogeneticTree:

	ParentItem = collections.namedtuple("ParentItem", ['name', 'distance'])

	def __init__(self, file, skipInit = False):

		if type(file) == tuple:
			print >> sys.stderr, "Creation of the phylogenetic tree ...",
			(self.items, self.root, self.officialName) = file
		else:
			print >> sys.stderr, "Loading phylogenetic tree %s ..." % file,
			self.officialName = {}
			self.items = self.newCommonNamesMapperInstance()

			# name and instance of file
			f = myFile.openFile(file, 'r')
			try:
				self.name = f.name
			except AttributeError:
				self.name = file

			f = myFile.firstLineBuffer(f)
			if (';' in f.firstLine) or ('(' in f.firstLine):
				self.__loadFromNewick__(' '.join(f).replace('\n','') + " ;")
			else:
				self.__loadFromMyFormat__(f)
			f.close()
			if not skipInit:
				self.reinitTree()
			else:
				print >> sys.stderr, "OK"


	def reinitTree(self):

		# analysing process of the tree
		def recInitialize(node):

			sys.stderr.write(".")

			self.dicLinks.setdefault(node, self.newCommonNamesMapperInstance())
			self.dicLinks.get(node).setdefault(node, [node])
			self.dicParents.setdefault(node, self.newCommonNamesMapperInstance())
			self.dicParents.get(node).setdefault(node, node)
			res = [node]

			if node in self.items:
				self.fileName.setdefault(node, str(node).replace(' ', '_').replace('/', '-'))
				s = []
				ld = []
				self.tmpA.append(node)
				for (f,l) in self.items.get(node):
					self.parent.setdefault(f, PhylogeneticTree.ParentItem(node,l))
					desc = recInitialize(f)
					s.extend(self.species.get(f))
					# we go back up
					res.extend(desc)
					ld.append(desc)
					for e in desc:
						self.dicParents.get(e).setdefault(node, node)
						self.dicParents.get(node).setdefault(e, node)
						self.dicLinks.get(e).setdefault(node, self.dicLinks.get(e).get(f) + [node])
						self.dicLinks.get(node).setdefault(e, [node] + self.dicLinks.get(f).get(e))
				self.species.setdefault(node, frozenset(s))
				# phylogenetic relationships
				for (s1,s2) in itertools.combinations(ld, 2):
					for e1 in s1:
						for e2 in s2:
							self.dicParents.get(e1).setdefault(e2, node)
							self.dicParents.get(e2).setdefault(e1, node)
							self.dicLinks.get(e1).setdefault(e2, self.dicLinks.get(e1).get(node) + self.dicLinks.get(node).get(e2)[1:])
							self.dicLinks.get(e2).setdefault(e1, self.dicLinks.get(e2).get(node) + self.dicLinks.get(node).get(e1)[1:])
			else:
				self.fileName.setdefault(node, str(node).replace(' ', '.'))
				self.species.setdefault(node, frozenset([node]))
				self.tmpS.append(node)

			self.allDescendants.setdefault(node, frozenset(res))

			return res

		print >> sys.stderr, "Data analysis ",

		# object initialisations
		self.parent = self.newCommonNamesMapperInstance()
		self.species = self.newCommonNamesMapperInstance()
		self.fileName = self.newCommonNamesMapperInstance()
		self.dicLinks = self.newCommonNamesMapperInstance()
		self.dicParents = self.newCommonNamesMapperInstance()
		self.allDescendants = self.newCommonNamesMapperInstance()
		self.tmpS = []
		self.tmpA = []

		# filling
		recInitialize(self.root)
		self.listSpecies = frozenset(self.species).difference(self.items)
		self.listAncestr = frozenset(self.items)

		# post-analysis
		self.allNames = self.tmpA + self.tmpS
		self.outgroupSpecies = self.newCommonNamesMapperInstance()
		self.indNames = self.newCommonNamesMapperInstance()
		for (i,e) in enumerate(self.allNames):
			self.outgroupSpecies.setdefault(e, self.listSpecies.difference(self.species.get(e)))
			self.indNames.setdefault(e, i)
			self.officialName[i] = e
		self.numItems = [[]] * len(self.allNames)
		for a in self.listAncestr:
			self.numItems[self.indNames.get(a)] = [(self.indNames.get(f),l) for (f,l) in self.items.get(a)]
		self.numParent = [None] * len(self.allNames)
		for (e,(par,l)) in self.parent.iteritems():
			self.numParent[self.indNames.get(e)] = (self.indNames.get(par),l)
		# official names / common names
		tmp = collections.defaultdict(list)
		for (curr,off) in self.officialName.iteritems():
			tmp[off].append(curr)
		self.commonNames = self.newCommonNamesMapperInstance()
		self.commonNames.update(tmp)
		self.dicGenes = {}
		self.dicGenomes = self.newCommonNamesMapperInstance()

		print >> sys.stderr, " OK"


	# return the name of the last common ancestor of several species
	def lastCommonAncestor(self, species):
		anc = species[0]
		for e in species[1:]:
			anc = self.dicParents[anc][e]
			if anc == self.root:
				return self.root
		return anc

	# assess if 'child' is really the child of 'parent'
	def isChildOf(self, child, parent):
		return self.dicParents[child][parent] == self.officialName[parent]
	
	#FIXME
	def compact(self, maxlength=1e-4):
		def do(node):
			if node not in self.items:
				return
			newitems = []
			for x in self.items[node]:
				do(x[0])
				if (x[1] >= maxlength) or (x[0] not in self.items):
					newitems.append(x)
				else:
					print >> sys.stderr, "Removing node", x[0]
					newitems.extend(self.items.pop(x[0]))
			self.items[node] = newitems
		do(self.root)
		self.reinitTree()

	# return the list of species to be compared
	# ancestors must be prefixed with '.' to be taken into account as genomes
	#   otherwise, that is their descendant species that are used 
	def getTargets(self, s):

		# are we exclusively on the common ancestors 
		allIntermediates = (s[-1] == "+")
		if allIntermediates:
			s = s[:-1]

		# list of extant species
		listSpecies = set()
		listAncestors = set()
		for x in s.split(','):
			if x[0] != '.':
				listSpecies.update(self.species[x])
			else:
				listSpecies.add(x[1:])
				listAncestors.add(x[1:])

		# list of ancestors
		for (e1,e2) in itertools.combinations(listSpecies, 2):
			if allIntermediates:
				listAncestors.update(self.dicLinks[e1][e2][1:-1])
			else:
				listAncestors.add(self.dicParents[e1][e2])
			# if e1 or e2 is already an ancestor

		return (listSpecies, listAncestors)

	# return the list of species pointed out by target 
	#   =esp
	#   +anc [children]
	#   /anc [outgroup]
	#   _** in order to remove instead of adding
	def getTargetsSpec(self, target):
		lesp = set()
		for x in target.split(","):
			if x.startswith("_"):
				if x.startswith("_/"):
					lesp.difference_update(self.outgroupSpecies[x[2:]])
				elif x.startswith("_="):
					lesp.discard(x[2:])
				else:
					lesp.difference_update(self.species[x[1:]])
			else:
				if x.startswith("/"):
					lesp.update(self.outgroupSpecies[x[1:]])
				elif x.startswith("="):
					lesp.add(x[1:])
				else:
					lesp.update(self.species[x])
		return lesp


	# return the list of ancestors pointed out by target
	#   =anc
	#   +anc [children]
	#   -anc [parents]
	#   /anc [parents+outgroups]
	#   _** in order to remove instead of adding
	def getTargetsAnc(self, target):
		lanc = set()
		for x in target.split(","):
			if x.startswith("_"):
				if x.startswith("_/"):
					lanc.difference_update(e for e in self.listAncestr if not self.isChildOf(e, x[2:]))
				elif x.startswith("_="):
					lanc.discard(x[2:])
				elif x.startswith("_\\"):
					lanc.difference_update(e for e in self.listAncestr if self.isChildOf(x[2:], e))
				else:
					lanc.difference_update(e for e in self.listAncestr if self.isChildOf(e, x[1:]))
			else:
				if x.startswith("/"):
					lanc.update(e for e in self.listAncestr if not self.isChildOf(e, x[1:]))
				elif x.startswith("="):
					lanc.add(x[1:])
				elif x.startswith("\\"):
					lanc.update(e for e in self.listAncestr if self.isChildOf(x[1:], e))
				else:
					lanc.update(e for e in self.listAncestr if self.isChildOf(e, x))
		return lanc


	# return the structure of the sub-tree in which are exclusively the chosen species
	def getSubTree(self, goodSpecies, rootAnc=None):
		goodAnc = set(self.dicParents[e1][e2] for (e1,e2) in itertools.combinations(goodSpecies, 2))
		newtree = collections.defaultdict(list)
		import myMaths
		def do(node):
			if node in goodAnc:
				for (x,_) in self.items[node]:
					newtree[node].extend(do(x))
				return [node]
			elif node in self.items:
				return myMaths.flatten([do(x) for (x,_) in self.items[node]])
			elif node in goodSpecies:
				return [node]
			else:
				return []
		root = do(self.root if rootAnc is None else rootAnc)
		assert len(root) == 1
		return (root[0], newtree)


	# calculate the values for the nodes of the tree
	#  - 'values' represents some defined values (nodes or leaves)
	#  - 'notdefined' is the value to be returned if no result
	#  - 'resultNode' indicates at which ancestral node we want to find the result
	def calcWeightedValue(self, values, notdefined, resultNode):

		import numpy

		n = len(self.allNames)
		# by default the results will be wil be "notdefined"
		matriceA = numpy.identity(n)
		matriceB = numpy.empty((n,))
		matriceB.fill(notdefined)

		# recursively build the matrix
		def recBuild(ianc, parent, length):
			anc = self.allNames[ianc]

			#FIXME recursive calls: children are selected OK (supposing that we are it ourself)
			items = [(e,p) for (e,p) in self.numItems[ianc] if recBuild(e, ianc, p)]
			# parent if available
			if parent != None:
				items.append( (parent,length) )

			if anc in values:
				# if we have a value, "x = val" is assigned
				matriceB[ianc] = values.get(anc)
				return True

			elif len(items) >= 2:
				# if it has enough neighbours, the equation is written
				s = 0.
				for (e,p) in items:
					p = 1./max(p,0.00001)
					matriceA[ianc][e] = p
					s += p
				matriceA[ianc][ianc] = -s
				matriceB[ianc] = 0
				return True
			else:
				return False

		# matrix building
		if len(values) == 0:
			return matriceB

		rootNode = self.indNames[self.lastCommonAncestor(values.keys())]
		recBuild(rootNode, None, 0)
		# equation solving
		res = numpy.linalg.solve(matriceA, matriceB)


		# search the neirest value around 'node' knowing that we come from 'origin'
		def searchBestValue(node, origin, tripLength):

			# Does node have a defined value
			val = res[node]
			if val != notdefined:
				return (tripLength,node,val)

			# Paths to go through
			test = self.numItems[node][:]
			par = self.numParent[node]
			if par != None:
				test.append(par)
			# Iteration
			best = (1e300,node,notdefined)
			for (e,l) in test:
				# In order to not loop
				if e == origin:
					continue
				# Break
				if best[0] < tripLength+l:
					continue
				tmp = searchBestValue(e, node, tripLength+l)
				if tmp[0] < best[0]:
					best = tmp
			return best


		if resultNode == None:
			return res
		resultNode = self.indNames[resultNode]
		best = searchBestValue(resultNode, None, 0)
		return (best[0],self.allNames[best[1]],best[2])



	# return a dict that uses internally official names of taxons
	#  but that can be used with common names
	def newCommonNamesMapperInstance(self):

		dsi = dict.__setitem__
		dgi = dict.__getitem__

		class commonNamesMapper(dict):

			def __getitem__(d, name):
				if name in self.officialName:
					return dgi(d, self.officialName[name])
				else:
					return dgi(d, name)

			def __setitem__(d, name, value):
				if name in self.officialName:
					dsi(d, self.officialName[name], value)
				else:
					dsi(d, name, value)

		return commonNamesMapper()

	# load all the species that come from an ancestor
	def loadAllSpeciesSince(self, ancestr, template, **args):
		if ancestr in self.species:
			l = self.species[ancestr]
		else:
			l = self.listSpecies
		self.loadSpeciesFromList(l, template, **args)

	# load all the outgroup species of an ancestor
	def loadAllSpeciesBefore(self, ancestr, template, **args):
		if ancestr in self.outgroupSpecies:
			l = self.outgroupSpecies[ancestr]
		else:
			l = self.listSpecies
		self.loadSpeciesFromList(l, template, **args)

	# load all the species of a list
	def loadSpeciesFromList(self, lst, template, storeGenomes = True):

		import myGenomes

		for esp in lst:
			esp = self.officialName[esp]
			g = myGenomes.Genome(template % self.fileName[esp])
			if storeGenomes:
				self.dicGenomes[esp] = g
			for (x,(c,i)) in g.dicGenes.iteritems():
				#self.dicGenes[x] = (esp, c, i)
				self.dicGenes[x] = GeneSpeciesPosition(esp, c, i)

	# return the number of genes in each species for a given family
	def findFamilyComposition(self, fam):

		score = self.newCommonNamesMapperInstance()
		score.update(dict.fromkeys(self.officialName,[]))

		for g in fam:
			if g in self.dicGenes:
				(e,_,_) = self.dicGenes[g]
				score[e].append(g)

		return score


	# topple over ("basculer" in french) the tree around a node
	def reroot(self, node, useOutgroups, newname=0):
		self.tmpItems = self.items.copy()
		if useOutgroups:
			anc = node
			while anc in self.parent:
				(par,l) = self.parent[anc]
				self.tmpItems[anc] = self.tmpItems[anc] + [(par,l)]
				self.tmpItems[par] = [x for x in self.tmpItems[par] if x[0] != anc]
				anc = par
		self.tmpItems[newname] = self.tmpItems[node]



	# load a phylogenetic tree into the phylTree format (with tabulations)
	def __loadFromMyFormat__(self, f):
		
		# loading process
		def loadFile():
			lignes = []
			for ligne in f:

				# the final "\n" is removed and we cut owing to the "\t"
				l = ligne.split('\t')

				# indentation level
				indent = 0
				for x in l:
					if x != '':
						break
					indent += 1

				# the triplet (indentation,names,age) is recorded
				names = [x.strip() for x in l[indent].split('|')]
				if len(l) > (indent+1):
					age = int(l[indent+1])
				else:
					age = 0
				lignes.append( (indent,names,age) )

			lignes.reverse()
			return lignes


		# lines analysis process
		def recLoad(indent):

			# we continue only if the next line is indented as expected
			if len(lignes) == 0  or lignes[-1][0] != indent:
				return None

			# the line is loaded
			currLine = lignes.pop()

			# all the children sub-trees are loaded, as long as it is possible
			children = []
			while True:
				tmp = recLoad(indent + 1)
				if tmp == None:
					break
				# record (name_of_child, evolution_time)
				children.append( (tmp, currLine[2]-self.ages.get(tmp)) )

			n = currLine[1][0]
			if len(children) == 0:
				if n[0] == ".":
					currLine[1][0] = n = n[1:]
					self.lstEsp6X.add(n)
				elif n[0] == "*":
					currLine[1][0] = n = n[1:]
					self.lstEsp2X.add(n)
				else:
					self.lstEspFull.add(n)

			# a single child, the node is returned
			if len(children) == 1:
				return children[0][0]
			# several children, they are recorded
			elif len(children) > 1:
				self.items.setdefault(n, children)

			# standard informations
			self.ages.setdefault(n, currLine[2])
			for s in currLine[1]:
				self.officialName[s] = n

			return n

		self.ages = self.newCommonNamesMapperInstance()
		self.lstEsp2X = set()
		self.lstEsp6X = set()
		self.lstEspFull = set()
		lignes = loadFile()
		self.root = recLoad(0)


	# tree in the Newick format (between parentheses)
	def __loadFromNewick__(self, s):

		#FIXME read the nb next characters of the tree
		def readStr(nb):
			ret = s[self.pos:self.pos+nb]
			self.pos += nb
			return ret

		# remove forbidden characters
		def keepWhile(car):
			x = self.pos
			while s[x] in car:
				x += 1
			return readStr(x-self.pos)

		# keep the authorised characters
		def keepUntil(car):
			x = self.pos
			while s[x] not in car:
				x += 1
			return readStr(x-self.pos)

		# read the tree in the form of text informations
		def readTree():
			keepWhile(' ')

			if s[self.pos] == '(':
				items = []
				# '(' the first time, then some ',' untill the final ')'
				while readStr(1) != ')':
					items.append( readTree() )
					keepWhile(' ')
				keepWhile(' ')
				# the result is the list of children + the name
				elt = (items, keepUntil("),:;[ "))
			else:
				# the result is a name
				elt = keepUntil("),:;[ ")

			keepWhile(' ')

			# possibly a non-null branch length
			if s[self.pos] == ':':
				readStr(1) # ":"
				length = float(keepWhile("0123456789.eE-"))
			else:
				length = 0

			keepWhile(' ')

			# possibly informations between brackets
			if s[self.pos] == '[':
				readStr(1) # "["
				info = keepUntil("]")
				info = dict( x.split("=") for x in info.split(":") if "=" in x )
				readStr(1) # "]"
			else:
				info = {}

			keepWhile(' ')

			return (elt,length,info)

		# fill the variables of a GenericTree
		def storeTree(data):
			(elt,length,info) = data

			if type(elt) == tuple:
				(children,nom) = elt
			else:
				children = []
				nom = elt

			if (nom == '') or (nom in self.officialName):
				nom = "NAME_%d" % self.pos
				self.pos += 1
			self.officialName[nom] = nom
			self.info[nom] = info

			if len(children) > 0:
				items = []
				for arbre in children:
					items.append( storeTree(arbre) )
				self.items.setdefault(nom, items)

			self.root = nom
			return (nom,length)

		self.pos = 0
		data = readTree()
		self.pos = 0
		self.info = {}
		storeTree(data)

