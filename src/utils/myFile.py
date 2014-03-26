# -*- coding: utf-8 -*-
# PhylDiag v1.01
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work
# under the terms of the GNU General Public License, version 2 or later

# Contient les fonctions de gestion des fichiers

import os

null = open(os.devnull, 'w')

###################################
# Gestion des fichiers tabulaires #
###################################
class myTSV:

	import collections
	csvProxy = collections.namedtuple("csvProxy", ['file','csvobject'])

	# Lecture en utilisant le module csv
	######################################
	@staticmethod
	def reader(fileName, **keywords):
		import csv
		f = openFile(fileName, 'r')
		return myTSV.csvProxy(f,csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n", **keywords))

	# Ecriture en utilisant le module csv
	#######################################
	@staticmethod
	def writer(fileName):
		import csv
		f = openFile(fileName, 'w')
		return myTSV.csvProxy(f,csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n"))

	# Renvoie la ligne preparee pour l'impression
	###############################################
	@staticmethod
	def printLine(line, delim = "\t", func = str):
		return delim.join(func(x) for x in line)


	# Lit un fichier tabulaire, en convertissant les colonnes separees de delim selon type_list
	#############################################################################################
	@staticmethod
	def readTabular(filename, type_list, delim = '\t'):

		import itertools
		f = openFile(filename, 'r')
		# Liste des types de chaque colonne
		new_type_list = []
		for x in type_list:
			if type(x) == type:
				new_type_list.append(x)
			else:
				new_type_list.extend([x[0]] * x[1])
		# Parcours du fichier
		for (i,line) in enumerate(f):
			current_line = line.replace('\n','').split(delim)
			assert len(current_line) == len(new_type_list), "Erreur nombre de colonne. Ligne:%d" % (i+1)
			yield tuple(t(x) for (x,t) in itertools.izip(current_line,new_type_list))
		f.close()

	# Permet de charger les dumps de MySQL (raboute les lignes tronquees)
	#######################################################################
	@staticmethod
	def MySQLFileLoader(f):
		tmp = ""
		for ligne in f:
			ligne = ligne.replace('\n', '')
			if ligne[-1] == '\\':
				# Signe que la ligne n'est pas terminee
				tmp = ligne[:-1]
			else:
				yield tmp + ligne
				tmp = ""
		assert (tmp == "")

	# Ecrit un fichier de dump MySQL (\N pour NULL)
	################################################
	@staticmethod
	def MySQLFileWriter(data):
		return myTSV.printLine(data).replace("None", "\N")

#########################################################################################################################
# Le but est de pouvoir acceder au fichier et lire la premiere ligne sans devoir le fermer pour le reouvrir juste apres #
#########################################################################################################################
class firstLineBuffer:
	def __init__(self, f):
		self.f = f
		try:
			self.firstLine = self.next()
		except StopIteration:
			self.firstLine = ""

	def __iter__(self):
		yield self.firstLine
		while True:
			yield self.next()
	
	def next(self):
		while True:
			l = self.f.next().replace('\n', '').replace('\r', '')
			# Suppression des lignes avec commentaire
			if (not l.startswith("#")) and (len(l) > 0):
				return l

	def close(self):
		return self.f.close()


####################
# Fichier existant #
####################
def hasAccess(s):
	return os.access(os.path.expanduser(s), os.R_OK)



####################################################################
# Cette fonction ouvre le fichier en le decompressant s'il le faut #
#   Retourne l'objet FILE et le nom complet du fichier             #
####################################################################
def openFile(nom, mode):

	# Fichier deja ouvert
	if type(nom) != str:
		return nom

	# Resource Web
	elif nom.startswith("http://") or nom.startswith("ftp://"):
		comm = "wget %s -O -"
		# Compression bzip2
		if nom.endswith(".bz2"):
			comm += " | bunzip2"
		# Compression gzip
		elif nom.endswith(".gz"):
			comm += " | gunzip"
		# Compression lzma
		elif nom.endswith(".lzma"):
			comm += " | unlzma"
		(stdin,f,stderr) = os.popen3( comm % nom )
		stdin.close()
		stderr.close()

	# Entree standard
	elif nom == "-":
		import sys
		return sys.stdin

	# Fichier sur le disque
	else:
		nom = os.path.expanduser(nom)
		if ("w" in mode) or ("a" in mode):
			# Cree le repertoire pour les sorties dans fichiers #
			try:
				os.makedirs(os.path.dirname(nom))
			except OSError:
				pass
		i = nom.find(".zip/")
		if (mode == "r") and (i >= 0):
			import zipfile
			import cStringIO
			f = zipfile.ZipFile(nom[:i+4], "r")
			f = cStringIO.StringIO(f.read(nom[i+5:]))
		# Compression bzip2
		elif nom.endswith(".bz2"):
			import bz2
			f = bz2.BZ2File(nom, mode)
		# Compression gzip
		elif nom.endswith(".gz"):
			import gzip
			f = gzip.GzipFile(nom, mode)
		# Compression lzma
		elif nom.endswith(".lzma") or nom.endswith(".xz"):
			import lzma
			f = lzma.LZMAFile(nom, mode)
		else:
			f = open(nom, mode)
	return f


