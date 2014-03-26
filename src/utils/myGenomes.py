# -*- coding: utf-8 -*-
# PhylDiag v1.01
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work
# under the terms of the GNU General Public License, version 2 or later

import sys
import itertools
import collections

import enum

import myFile
import myTools


Gene = collections.namedtuple("Gene", ['chromosome', 'beginning', 'end', 'strand', 'names'])
GenePosition = collections.namedtuple("GenePosition", ['chromosome', 'index'])


###########################################################
# On convertit les noms de chromosomes en int si possible #
###########################################################
def commonChrName(x):
	try:
		return int(x)
	except TypeError:
		return None
	except ValueError:
		return intern(x)


############################################
# Le type de contig en fonction de son nom #
############################################
ContigType = enum.Enum('None', 'Random', 'Scaffold', 'Chromosome')
def contigType(chrom):
	if chrom in [None, "Un_random", "UNKN", "Un"]:
		return ContigType.None
	try:
		x = int(chrom)
		if x < 100:
			return ContigType.Chromosome
		else:
			return ContigType.Scaffold
	except:
		c = chrom.lower()
		if "rand" in c:
			return ContigType.Random
		keys = ["cont", "scaff", "ultra", "reftig", "_", "un", "mt", "gl"]
		for x in keys:
			if x in c:
				return ContigType.Scaffold
		if (chrom in ["U", "E64", "2-micron"]) or chrom.endswith("Het"):
			return ContigType.Scaffold
		else:
			return ContigType.Chromosome


class myFASTA:

	# Charge un fichier FASTA une sequence a la fois
	#################################################
	@staticmethod
	def iterFile(name):
		f = myFile.openFile(name, "r")
		name = None
		for ligne in f:
			ligne = ligne.replace('\n', '').strip()
			# Les chevrons indiquent le debut d'une nouvelle sequence
			if ligne.startswith('>'):
				if name != None:
					yield (name, "".join(tmp))
				name = ligne[1:].strip()
				tmp = []
			# Les lignes doivent etre concatenees
			elif name != None:
				tmp.append(ligne.upper())
		if name != None:
			yield (name, "".join(tmp))
		f.close()

	# Charge un fichier FASTA en entier
	####################################
	@staticmethod
	def loadFile(name):
		return dict(myFASTA.iterFile(name))

	# Ecrit une sequence au format FASTA
	#####################################
	@staticmethod
	def printSeq(f, name, seq, length=60):
		print >> f, ">" + name
		while len(seq) != 0:
			print >> f, seq[:length]
			seq = seq[length:]
		print >> f



# Comptage de nucleotides
##########################
def getMonoNuc(seq, x1, x2, sel):
	n = 0
	gc = 0
	for x in xrange(x1, x2+1):
		if seq[x] == "N":
			continue
		n += 1
		if seq[x] in sel:
			gc += 1
	return (n,gc)

# Comptage de dinucleotides
############################
def getDiNuc(seq, x1, x2, sel):
	n = 0
	gc = 0
	for x in xrange(x1, x2):
		s = seq[x:x+2]
		if "N" in s:
			continue
		n += 1
		if s in sel:
			gc += 1
	return (n,gc)


codon2aa = {
     'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
     'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
     'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
     'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
     'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*' }

aa2codon = collections.defaultdict(list)
for (c,p) in codon2aa.iteritems():
	aa2codon[p].append(c)



##########################################
# Classe generale pour stocker un genome #
##########################################
class Genome:

	# Constructeur
	###############
	def __init__(self, fichier, **kwargs):

		if isinstance(fichier, str):
			print >> sys.stderr, "Chargement du genome de", fichier, "...",

			f = myFile.firstLineBuffer(myFile.openFile(fichier, 'r'))

			# la liste des genes par chromosome
			self.lstGenes = collections.defaultdict(list)

			# Le choix de la fonction de chargement
			c = f.firstLine.split("\t")
			if f.firstLine.startswith(">") or f.firstLine.endswith("$"):
				# Format GRIMM Synteny
				#######################
				if f.firstLine.startswith(">"):
					self.name = f.firstLine[1:].strip()
				chrom = 1
				for l in f:
					l = l.strip()
					if not l.endswith("$"):
						continue
					for (i,x) in enumerate(l.replace("$","").split()):
						strand = -1 if x.startswith("-") else 1
						self.addGene([x[1:] if x[0] in "-+" else x], chrom, i, i+1, strand)
					chrom += 1
				print >> sys.stderr, "(GRIMM)",

			elif len(c) == 1:

				# Genes ancestraux: "NAMES"
				#############################
				for (i,l) in enumerate(f):
					self.lstGenes[None].append( Gene(None, i, i+1, 0, tuple(intern(x) for x in l.split())) )
				print >> sys.stderr, "(genes ancestraux)",

			elif (len(c) == 2) and not set(c[1]).issubset("01-"):

				# Genome ancestral: "CHR NAMES"
				#################################
				lastC = None
				for (i,l) in enumerate(f):
					c = l.split("\t")
					if lastC != c[0]:
						lastC = c[0]
						dec = i
					self.addGene(c[1].split(), c[0], i-dec, i-dec+1, 0)
				print >> sys.stderr, "(genome ancestral: chrom+noms)",

			elif (len(c) >= 5) and (" " not in c[3]) and (len(c[4]) > 0):

				# Ensembl: "CHR BEG END STRAND NAMES"
				#################################################
				for l in f:
					c = l.replace('\n', '').split('\t')
					self.addGene(c[4].split(), c[0], int(c[1]), int(c[2]), int(c[3]) if c[3]!='None' else None)
				print >> sys.stderr, "(Ensembl)",

			elif (len(c) == 4) and int(c[1]) < 2:

				# Genome ancestral: "CHR STRAND LST-INDEX LST-STRANDS"
				########################################################
				
				if 'ancGenes' in kwargs:
					ancGenes = kwargs["ancGenes"].lstGenes[None]

				lastC = None
				for l in f:
					c = l.split("\t")
					if lastC != c[0]:
						lastC = c[0]
						pos = 0
						currC = commonChrName(c[0])
					data = zip([int(x) for x in c[2].split()], [int(x) for x in c[3].split()])
					if int(c[1]) < 0:
						data = [(i,-s) for (i,s) in data.__reversed__()]
					for (index,strand) in data:
						if 'ancGenes' in kwargs:
							self.lstGenes[currC].append( Gene(currC, pos, pos+1, strand, ancGenes[index].names) )
						else:
							self.lstGenes[currC].append( Gene(currC, pos, pos+1, strand, (index,)) )
						pos += 1
				print >> sys.stderr, "(genome ancestral: chrom+diags)",

			else:
				if len(c) == 2:
					(ili,ils) = (0,1)
				else:
					assert len(c) >= 4
					(ili,ils) = (2,3)
					self.ancName = c[0]

				if 'ancGenes' in kwargs:
					ancGenes = kwargs["ancGenes"].lstGenes[None]

				# Genome ancestral: "LST-INDEX LST-STRANDS"
				#############################################
				for (i,l) in enumerate(f):
					c = l.split("\t")
					chrom = i+1
					lchrom = self.lstGenes[chrom]
					for (pos,(index,strand)) in enumerate(itertools.izip(c[ili].split(), c[ils].split())):
						if 'ancGenes' in kwargs:
							lchrom.append( Gene(chrom, pos, pos+1, int(strand), ancGenes[int(index)].names) )
						else:
							lchrom.append( Gene(chrom, pos, pos+1, int(strand), (int(index),) ) )
				print >> sys.stderr, "(genome ancestral: diags)",

			f.close()
			self.name = fichier

		else:
			genomeBase = fichier
			print >> sys.stderr, "Filtrage de", genomeBase.name, "...",
			filterIn = set(kwargs["filterIn"]) if "filterIn" in kwargs else None
			filterOut = set(kwargs["filterOut"]) if "filterOut" in kwargs else None

			def filt(gene):
				if filterIn is not None:
					return any(s in filterIn for s in gene.names)
				if filterOut is not None:
					return all(s not in filterOut for s in gene.names)
				return True

			self.lstGenes = {}
			for (chrom,l) in genomeBase.lstGenes.iteritems():
				l = [gene for gene in l if filt(gene)]
				if len(l) > 0:
					self.lstGenes[chrom] = l
			self.name = "Filter from " + genomeBase.name
			print >> sys.stderr, "%d genes -> %d genes" % (sum(len(x) for x in genomeBase.lstGenes.itervalues()), sum(len(x) for x in self.lstGenes.itervalues())),

		self.init(**kwargs)
		print >> sys.stderr, "OK"


	# Initialise le dictionnaire et la liste des chromosomes
	##########################################################
	def init(self, **kwargs):

		self.dicGenes = collections.defaultdict(list)
		self.chrList = collections.defaultdict(list)
		self.chrSet = collections.defaultdict(set)
		withDict = kwargs.get("withDict", True)

		for chrom in self.lstGenes:

			# Associer un nom de gene a sa position sur le genome
			self.lstGenes[chrom].sort()
			if withDict:
				for (i,gene) in enumerate(self.lstGenes[chrom]):
					for s in gene.names:
						self.dicGenes[s].append(GenePosition(chrom,i))

			# Classification en chromosomes/contigs/scaffolds
			t = contigType(chrom)
			self.chrList[t].append(chrom)
			self.chrSet[t].add(chrom)

		# Trie les noms des chromosomes
		for t in self.chrList:
			self.chrList[t].sort()


	# Rajoute un gene au genome
	############################
	def addGene(self, names, chromosome, beg, end, strand):
		#print >> sys.stderr, names
		assert 0 <= beg <= end
		assert strand in [-1,0,1,None]

		chromosome = commonChrName(chromosome)
		self.lstGenes[chromosome].append( Gene(chromosome, beg, end, strand, tuple(intern(s) for s in names)) )


	# Renvoie les genes presents sur le chromosome donne entre beg et end paires de bases
	#####################################################################################
	def getGenesAt(self, chr, beg, end, onlyInside=False):
		if chr not in self.lstGenes:
			return
		lst = self.lstGenes[chr]

		# Recherche dichotomique d'un gene d'intersection non vide avec la fenetre
		def dichotFind(a, b):
			i = (a+b)/2
			if a == i: # si b == a+1
				return a
			if lst[i].end < beg:
				return dichotFind(i,b)
			elif lst[i].beginning > end:
				return dichotFind(a,i)
			else:
				return i
		index = dichotFind(0, len(lst)-1)

		res = []
		for i in xrange(index-1, -1, -1):
			g = lst[i]
			if g.end < beg:
				break
			#if g.beginning <= end: # Toujours vrai
			assert g.beginning <= end
			if (not onlyInside) or ((beg <= g.beginning) and (g.end <= end)):
				res.append(g)
		for g in reversed(res):
			yield g

		for i in xrange(index, len(lst)):
			g = lst[i]
			if g.beginning > end:
				break
			#if g.end >= beg: # Toujours vrai
			assert g.end >= beg
			if (not onlyInside) or ((beg <= g.beginning) and (g.end <= end)):
				yield g



	# Renvoie les genes presents aux alentours d'un gene donne (fenetre de + ou - l en nombre de genes autour du gene specifie)
	############################################################################################################################
	def getGenesNear(self, chr, index, l):
		if chr not in self.lstGenes:
			return []

		return self.lstGenes[chr][max(0, index-l):index+l+1]


	# Renvoie tous les genes
	#########################
	def __iter__(self):
		for t in self.lstGenes.itervalues():
			for g in t:
				yield g


	# Cherche les positions d'un gene donne par ses noms
	####################################################
	def getPositions(self, names):
		return set( (gene for s in names if s in self.dicGenes for gene in self.dicGenes[s]) )

	# Cherche une position du gene donne par ses noms
	####################################################
	def getPosition(self, names):
		return set( (self.dicGenes[s][0] for s in names if s in self.dicGenes) )


	# Renvoie les autres noms d'un gene (en supposant que celui ci ne se trouve qu'a un endroit)
	####################################
	def getOtherNames(self, name):
		if name not in self.dicGenes:
			return []
		(c,i) = self.dicGenes[name][0]
		return [x for x in self.lstGenes[c][i].names if x != name]


	# Fabrique la liste des orthologues entre les deux genomes
	###########################################################
	def buildOrthosTable(self, chr1, genome2, chr2, includeGaps, genesAnc):

		chr2 = set(chr2)
		# Tous les orthologues entre pour les chromosomes OK du genome 1
		res = {}
		for c1 in chr1:
			res[c1] = []
			for (i1,g1) in enumerate(self.lstGenes[c1]):
				tmp = genome2.getPosition(g1.names)
				if genesAnc != None:
					for (c,i) in genesAnc.getPosition(g1.names):
						tmp.update(genome2.getPosition(genesAnc.lstGenes[c][i].names))

				# On ne garde que les chromsomes OK du genome 2
				tmp = [(c,i) for (c,i) in tmp if c in chr2]
				# +/- includeGaps
				if includeGaps or (len(tmp) > 0):
					res[c1].append( (i1,tmp) )
		return res


	# Imprime le genome (au format Ensembl)
	########################################
	def printEnsembl(self, f):
		for chrom in self.lstGenes:
			for gene in self.lstGenes[chrom]:
				print >> f, myFile.myTSV.printLine([chrom, gene.beginning, gene.end, gene.strand, " ".join(gene.names)])

