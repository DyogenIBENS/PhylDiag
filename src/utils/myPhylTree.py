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


##########################################
# Classe mere de tous les types d'arbres #
##########################################
class PhylogeneticTree:

	ParentItem = collections.namedtuple("ParentItem", ['name', 'distance'])

	def __init__(self, fichier, skipInit = False):

		if type(fichier) == tuple:
			print >> sys.stderr, "Creation de l'arbre phylogenetique ...",
			(self.items, self.root, self.officialName) = fichier
		else:
			print >> sys.stderr, "Chargement de l'arbre phylogenique de %s ..." % fichier,
			self.officialName = {}
			self.items = self.newCommonNamesMapperInstance()

			# le nom et l'instance de file
			f = myFile.openFile(fichier, 'r')
			try:
				self.nom = f.name
			except AttributeError:
				self.nom = fichier

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

		# La procedure d'analyse de l'arbre
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
					# On remonte a node
					res.extend(desc)
					ld.append(desc)
					for e in desc:
						self.dicParents.get(e).setdefault(node, node)
						self.dicParents.get(node).setdefault(e, node)
						self.dicLinks.get(e).setdefault(node, self.dicLinks.get(e).get(f) + [node])
						self.dicLinks.get(node).setdefault(e, [node] + self.dicLinks.get(f).get(e))
				self.species.setdefault(node, frozenset(s))
				# Liens de parente
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

		print >> sys.stderr, "Analyse des donnees ",

		# Initialisation des structures
		self.parent = self.newCommonNamesMapperInstance()
		self.species = self.newCommonNamesMapperInstance()
		self.fileName = self.newCommonNamesMapperInstance()
		self.dicLinks = self.newCommonNamesMapperInstance()
		self.dicParents = self.newCommonNamesMapperInstance()
		self.allDescendants = self.newCommonNamesMapperInstance()
		self.tmpS = []
		self.tmpA = []

		# Remplissage
		recInitialize(self.root)
		self.listSpecies = frozenset(self.species).difference(self.items)
		self.listAncestr = frozenset(self.items)

		# Structures post-analyse
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
		# Les noms officiels / courants
		tmp = collections.defaultdict(list)
		for (curr,off) in self.officialName.iteritems():
			tmp[off].append(curr)
		self.commonNames = self.newCommonNamesMapperInstance()
		self.commonNames.update(tmp)
		self.dicGenes = {}
		self.dicGenomes = self.newCommonNamesMapperInstance()

		print >> sys.stderr, " OK"


	# Renvoie le nom de l'ancetre commun de plusieurs especes
	##########################################################
	def lastCommonAncestor(self, species):
		anc = species[0]
		for e in species[1:]:
			anc = self.dicParents[anc][e]
			if anc == self.root:
				return self.root
		return anc

	# Teste si le fils est bien un fils de son pere
	################################################
	def isChildOf(self, fils, pere):
		return self.dicParents[fils][pere] == self.officialName[pere]

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

	#
	# Forme la liste des especes a comparer
	# Les ancetres doivent etre prefixes d'un '.' pour etre pris en tant que genomes
	#   sinon, ce sont leurs especes descendantes qui sont utilisees
	##################################################################################
	def getTargets(self, s):

		# Est-ce qu'on est uniquement sur les ancetres communs
		allIntermediates = (s[-1] == "+")
		if allIntermediates:
			s = s[:-1]

		# La liste des especes modernes
		listSpecies = set()
		listAncestors = set()
		for x in s.split(','):
			if x[0] != '.':
				listSpecies.update(self.species[x])
			else:
				listSpecies.add(x[1:])
				listAncestors.add(x[1:])

		# La liste des ancetres
		for (e1,e2) in itertools.combinations(listSpecies, 2):
			if allIntermediates:
				listAncestors.update(self.dicLinks[e1][e2][1:-1])
			else:
				listAncestors.add(self.dicParents[e1][e2])
			# si e1 ou e2 est deja un ancetre

		return (listSpecies, listAncestors)

	#
	# Renvoie la liste des especes designees par target
	#   =esp
	#   +anc [enfants]
	#   /anc [outgroup]
	#   _** pour enlever au lieu de rajouter
	#######################################################
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


	#
	# Renvoie la liste des ancetres designes par target
	#   =anc
	#   +anc [enfants]
	#   -anc [parents]
	#   /anc [parents+outgroups]
	#   _** pour enlever au lieu de rajouter
	##################################################################################
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


	# Renvoie la structure du sous-arbre dans lequel se trouve uniquement les especes choisies
	###########################################################################################
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


	# Calcule les valeurs sur les noeuds de l'arbre
	#  - values represente des valeurs definies (noeuds ou feuilles)
	#  - notdefined est la valeur a renvoyer si pas de resultat
	#  - resultNode indique de quel ancetre on veut le resultat
	#########################################################################
	def calcWeightedValue(self, values, notdefined, resultNode):

		import numpy

		n = len(self.allNames)
		# Par defaut, les resultats seront "notdefined"
		matriceA = numpy.identity(n)
		matriceB = numpy.empty((n,))
		matriceB.fill(notdefined)

		# Construit recursivement la matrice
		def recBuild(ianc, father, length):
			anc = self.allNames[ianc]

			# Appels recursifs: on selectionne les fils OK (en supposant qu'on le soit nous meme)
			items = [(e,p) for (e,p) in self.numItems[ianc] if recBuild(e, ianc, p)]
			# Le pere si disponible
			if father != None:
				items.append( (father,length) )

			if anc in values:
				# Si on a une valeur, on met "x = val"
				matriceB[ianc] = values.get(anc)
				return True

			elif len(items) >= 2:
				# S'il a suffisament de voisins, on ecrit l'equation
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

		# Construction de la matrice
		if len(values) == 0:
			return matriceB

		rootNode = self.indNames[self.lastCommonAncestor(values.keys())]
		recBuild(rootNode, None, 0)
		# Resolution de l'equation
		res = numpy.linalg.solve(matriceA, matriceB)


		# Cherche la valeur la plus proche de node sachant que notre parcours nous fait venir de origin
		def searchBestValue(node, origin, tripLength):

			# Est-ce que node a une valeur definie
			val = res[node]
			if val != notdefined:
				return (tripLength,node,val)

			# Les chemins a parcourir
			test = self.numItems[node][:]
			par = self.numParent[node]
			if par != None:
				test.append(par)
			# Iteration
			best = (1e300,node,notdefined)
			for (e,l) in test:
				# Pour ne pas boucler
				if e == origin:
					continue
				# Coupure
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



	# Renvoie un dictionnaire qui utilise en interne les noms officiels des taxons
	#  mais que l'on peut utiliser avec les noms communs
	# #############################################################################
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

	# Charge toutes les especes qui descendent d'un ancetre
	def loadAllSpeciesSince(self, ancestr, template, **args):
		if ancestr in self.species:
			l = self.species[ancestr]
		else:
			l = self.listSpecies
		self.loadSpeciesFromList(l, template, **args)

	# Charge toutes les especes outgroup d'un ancetre
	def loadAllSpeciesBefore(self, ancestr, template, **args):
		if ancestr in self.outgroupSpecies:
			l = self.outgroupSpecies[ancestr]
		else:
			l = self.listSpecies
		self.loadSpeciesFromList(l, template, **args)

	# Charge toutes les especes d'une liste
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

	# Renvoie le nombre de genes dans chaque espece pour une famille donnee
	def findFamilyComposition(self, fam):

		score = self.newCommonNamesMapperInstance()
		score.update(dict.fromkeys(self.officialName,[]))

		for g in fam:
			if g in self.dicGenes:
				(e,_,_) = self.dicGenes[g]
				score[e].append(g)

		return score


	#
	# Bascule l'arbre autour d'un noeud
	#
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



	# Charge un arbre phylogenetique dans mon format a moi
	#######################################################
	def __loadFromMyFormat__(self, f):

		# Procedure de chargement du fichier
		def loadFile():
			lignes = []
			for ligne in f:

				# On enleve le \n final et on coupe suivant les \t
				l = ligne.split('\t')

				# Le niveau d'indentation
				indent = 0
				for x in l:
					if x != '':
						break
					indent += 1

				# On stocke le triplet (indentation,noms,age)
				noms = [x.strip() for x in l[indent].split('|')]
				if len(l) > (indent+1):
					age = int(l[indent+1])
				else:
					age = 0
				lignes.append( (indent,noms,age) )

			lignes.reverse()
			return lignes


		# La procedure d'analyse des lignes du fichier
		def recLoad(indent):

			# On ne continue que si la ligne suivante est indentee comme prevu
			if len(lignes) == 0  or lignes[-1][0] != indent:
				return None

			# On charge la ligne
			currLine = lignes.pop()

			# On charge tous les sous arbres-fils tant que possible
			fils = []
			while True:
				tmp = recLoad(indent + 1)
				if tmp == None:
					break
				# On stocke (nom_du_fils, temps_d_evolution)
				fils.append( (tmp, currLine[2]-self.ages.get(tmp)) )

			n = currLine[1][0]
			if len(fils) == 0:
				if n[0] == ".":
					currLine[1][0] = n = n[1:]
					self.lstEsp6X.add(n)
				elif n[0] == "*":
					currLine[1][0] = n = n[1:]
					self.lstEsp2X.add(n)
				else:
					self.lstEspFull.add(n)

			# Un seul fils, on remonte le noeud
			if len(fils) == 1:
				return fils[0][0]
			# Plusieurs fils, on les enregistre
			elif len(fils) > 1:
				self.items.setdefault(n, fils)

			# Info standard
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


	# Arbre au format Newick (parentheses)
	#######################################
	def __loadFromNewick__(self, s):

		# Lit les nb prochains caracteres de l'arbre
		def readStr(nb):
			ret = s[self.pos:self.pos+nb]
			self.pos += nb
			return ret

		# Enleve les caracteres interdits
		def keepWhile(car):
			x = self.pos
			while s[x] in car:
				x += 1
			return readStr(x-self.pos)

		# Garde les caracteres autorises
		def keepUntil(car):
			x = self.pos
			while s[x] not in car:
				x += 1
			return readStr(x-self.pos)

		# Lit l'arbre sous forme de texte
		def readTree():
			keepWhile(' ')

			if s[self.pos] == '(':
				items = []
				# '(' la premiere fois, puis des ',' jusqu'au ')' final
				while readStr(1) != ')':
					items.append( readTree() )
					keepWhile(' ')
				keepWhile(' ')
				# Le resultat est la liste des fils + le nom
				elt = (items, keepUntil("),:;[ "))
			else:
				# Le resultat est un nom
				elt = keepUntil("),:;[ ")

			keepWhile(' ')

			# Eventuellement une longueur de branche non nulle
			if s[self.pos] == ':':
				readStr(1) # ":"
				length = float(keepWhile("0123456789.eE-"))
			else:
				length = 0

			keepWhile(' ')

			# Eventuellement des infos entre crochets
			if s[self.pos] == '[':
				readStr(1) # "["
				info = keepUntil("]")
				info = dict( x.split("=") for x in info.split(":") if "=" in x )
				readStr(1) # "]"
			else:
				info = {}

			keepWhile(' ')

			return (elt,length,info)

		# Remplit les variables d'un GenericTree
		def storeTree(data):
			(elt,length,info) = data

			if type(elt) == tuple:
				(fils,nom) = elt
			else:
				fils = []
				nom = elt

			if (nom == '') or (nom in self.officialName):
				nom = "NAME_%d" % self.pos
				self.pos += 1
			self.officialName[nom] = nom
			self.info[nom] = info

			if len(fils) > 0:
				items = []
				for arbre in fils:
					items.append( storeTree(arbre) )
				self.items.setdefault(nom, items)

			self.root = nom
			return (nom,length)

		self.pos = 0
		data = readTree()
		self.pos = 0
		self.info = {}
		storeTree(data)

