#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.01
# needs python 2.7 at least
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work
# under the terms of the GNU General Public License, version 2 or later

###################################################################
# PhylDiag core algorithm and tools to manage synteny blocks #
###################################################################

# This code is based on the publication ``PhylDiag : identifying complex cases of conserved synteny that include tandem duplications''
# authors :
# Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
# Matthieu MUFFATO (EBI, Cambridge, United Kingdom, muffato@ebi.ac.uk)
# Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)

# This code uses personal libraries : myTools and myGenomes

###################################
# Vocabulary (see the publication)
###################################
# g : gene
# sb : synteny block
# tb : tandem block
# s : may refer to 'strand', the transcriptional orientation of a gene or a tb
# hp : homology pack
# sign : may refer to a sign of a homology pack
# MH : matrix of homologies 
# MHP: matrix of homology packs
# diagonal : synonym of a sb (a sb is roughly  a daigonal in the MHP)
# diagonal type : either 'slash' ('/', bottom-left to top-right) or 'backslash' ('\', top-left to bottom-right)
# consistent diagonal : either a slash diagonal with hps signs = (+1 or None) or a backslash diagonal with hps signs = (-1 or None)

import os
import sys
import math
import copy
import itertools
import collections

import enum

import utils.myTools
import utils.myGenomes
import utils.myProbas

# For parallel computation
from multiprocessing import Queue, Process, Manager, Lock
import multiprocessing

FilterType = enum.Enum('None', 'InCommonAncestor', 'InBothSpecies')

# The frame of all the distances used during the merging process
def frameDistance(f,(x0,y0),(x1,y1), diagType):
	if diagType == '/':
		# Only consider the top right part of the distance matrix, (x1>=x0 and y1>=y0)
		if (x1-x0)>=0 and (y1-y0) >= 0: # in i-adhore "...>0", here we allow '=' to take care of close dispered paralogies that can be understood as tandem duplicates
			res = f((x0,y0),(x1,y1))
		else:
			res = sys.maxint
	elif diagType == '\\':
		# Only consider the bottom right part of the distance matrix, (x1>=x0 and y1<=y0)
		if (x1-x0) >= 0 and (y1-y0) <= 0: # in i-adhore "...<0", here we allow '=', idem than before
			res = f((x0,y0),(x1,y1))
		else:
			res = sys.maxint
	else:
		# diagType == None, diag is composed of a single gene
		res = f((x0,y0),(x1,y1))
	return res

# Diagonal Pseudo Distance
def DPD((x0,y0),(x1,y1), diagType):
	f=lambda (x0,y0),(x1,y1) : 2*max(abs(x1-x0),abs(y1-y0))-min(abs(x1-x0),abs(y1-y0))
	return frameDistance(f,(x0,y0),(x1,y1), diagType)

# Chebyshev Distance 
def CD((x0,y0),(x1,y1), diagType):
	f=lambda (x0,y0),(x1,y1) : max(abs(x1-x0), abs(y1-y0)) 
	return frameDistance(f,(x0,y0),(x1,y1), diagType)

# Manhattan Distance
def MD((x0,y0),(x1,y1), diagType):
	f=lambda (x0,y0),(x1,y1) : abs(x1-x0) + abs(y1-y0)
	return frameDistance(f,(x0,y0),(x1,y1), diagType)

# Euclidean Distance
def ED((x0,y0),(x1,y1), diagType):
	f=lambda (x0,y0),(x1,y1) : round(math.sqrt((x1-x0)**2 + (y1-y0)**2))
	return frameDistance(f,(x0,y0),(x1,y1), diagType)

#
# Generator managing the queue of diagonals for the merging process 
####################################################################
class queueWithBackup:
	# gen is a generator
	def __init__(self, gen):
		self.gen = gen
		self.backup = collections.deque()
		self.todofirst = []
	
	def __iter__(self):
		return self
	
	# The next returned value comes from either the buffer or from the main generator
	def next(self):
		if len(self.todofirst) > 0:
			return self.todofirst.pop()
		return self.gen.next() # Returns the next value of the generator
	
	# The reinserted element is put on hold, the last reinserted element will be on top of the pile during the next rewind
	def putBack(self, x):
		self.backup.appendleft(x)
	
	# Recorded elements are put in a prioritary buffer
	def rewind(self):
		self.todofirst.extend(self.backup)
		self.backup = collections.deque()

#
# Merge diagonals if they are separated by a gap less long than 'gapMax' relatively to a distance metric 
# inputs:
# listOfDiags : a list containing elements as (l1, l2, la)
#	l1 = [ (i12,s12), (i12,s12), ....] with i1x < i1y for all x<y
#		i1 : gene index in the genome i1 
#	l2 : [..., (i2,s2), ...] 
#		i2 : gene index in the genome i2 
#	la : [..., j1,...] 
#		j1 is the ancGene number (corresponds to the line index of the gene)  
# gapMax : the maximum allowed gap between merged diagonals
# distanceMetric : the metric used for the gap calculation
# outputs: 
#	listOfSortedAndMergedDiagonals : same structure than listOfDiags
########################################################################################
# TODO : optimize the search of diags by avoiding considering diags that are on the left of diagA
@utils.myTools.verbose
def mergeSbs(listOfDiags, gapMax, distanceMetric = 'DPD', verbose = True): 
	assert gapMax>=0

	if distanceMetric == 'DPD':
		print >> sys.stderr, "use Diagonal Pseudo Distance to merge diagonals with a gap up to %s elements" % gapMax 
		distance = DPD
	elif distanceMetric == 'CD':
		print >> sys.stderr, "use Chebyshev Distance to merge diagonals with a gap up to %s elements" % gapMax 
		distance = CD
	elif distanceMetric == 'MD':
		print >> sys.stderr, "use Manhattan Distance to merge diagonals with a gap up to %s elements" % gapMax 
		distance = MD
	elif distanceMetric == 'ED':
		print >> sys.stderr, "use Euclidean Distance to merge diagonals with a gap up to %s elements" % gapMax 
		distance = ED
	else:
		raise ErrorValue('Must use a distance either DPD (Diagonal PSeudo Distance) or MD (Manhattan Distance), Euclidean Distance (ED) or Chebyshev Distance (CD)') 

	print >> sys.stderr, "Number of Diags before DiagMerger = ", len(listOfDiags)
	diagGen = []
	listOfFinishedDiags = []
	nbFusion = 0
	listOfDiags = [(s,l1,l2,la,(l1[0][0],l2[0][0]), (l1[-1][0],l2[-1][0])) for (s,l1,l2,la) in listOfDiags]
	for currGap in range(0,gapMax+1):
		nbFusionCurrGap = 0
		# Sort diagonals by increasing index on the genome A (i.e. by increasing x coordinate of the leftmost hp of the diagonal)
		# This is a key part of the algorithm ! During the whole process we need to keep diagonals sorted by increasing index on the genome A
		if currGap > 0: #After initialisation
			listOfDiags = sorted(list(diagGen) + listOfFinishedDiags, key=lambda x:x[1][0])
			listOfFinishedDiags=[]
		# assert all([l1_1[0][0] <= l1_2[0][0] for ((_,l1_1,_,_,_,_),(_,l1_2,_,_,_,_)) in utils.myTools.myIterator.slidingTuple(listOfDiags)])
		# assert all([l1[0][0] == min(l1)[0] for (_,l1,_,_,_,_) in listOfDiags])
		# assert all([l1[-1][0] == max(l1)[0] for (_,l1,_,_,_,_) in listOfDiags])
		diagGen = queueWithBackup( d for d in listOfDiags )
		for diagA in diagGen:
			(dTa,la1,la2,laa,startA,endA) = diagA
			### DEBUG example
			#if 426 == startA[0]+1 and 543 == startA[1]+1:
			# 	print >> sys.stderr, "Is diag A : [diagType = %s, start = (%s,%s), end =(%s,%s)] fusionable with" % (dTa,startA[0]+1,startA[1]+1,endA[0]+1,endA[1]+1)
			### DEBUG example
			continueToSearch = True
			fusionableDiags = []
			impossibleToMergeDiags = []
			for diagB in diagGen:
				if not continueToSearch:
					impossibleToMergeDiags.append(diagB)
					break
				(dTb,_,_,_,startB,endB) = diagB
				### DEBUG example
				#if 427 == startB[0]+1 and 541 == startB[1]+1:
				#	print >> sys.stderr, "Is diag B : [diagType = %s, start = (%s,%s), end =(%s,%s)] fusionable with" % (dTb,startB[0]+1,startB[1]+1,endB[0]+1,endB[1]+1) 
				#	print >> sys.stderr, "endA = (%s,%s)" % (endA[0]+1,endA[1]+1), " ?"
				### DEBUG example
				# Thanks to the sort at the beginning, we known that the starting hp of diagB is on the right of the starting hp of diagA
				#TODO : change diagGen to put diags that are on the left of diagA in a buffer and avoid considering them 
				if startB[0] < endA[0]: # in i-adhore startB[0] <= endA[0]
					impossibleToMergeDiags.append(diagB)
					continue
				elif startB[0] <= endA[0] + currGap+1: 
					# endA[0] < startB[0] <= endA[0] + currGap
					# Check if diagTypes are compatible
					if dTa == dTb or dTa == None or dTb == None:
						# take the known diagType if it is known in at least one of the 2 diags
						dT = dTa if dTa != None else dTb
						# Check the distance
						if distance(endA,startB,dT) == currGap+1 :
								fusionableDiags.append(diagB)
						else:
							impossibleToMergeDiags.append(diagB)
							continue
					else:
						impossibleToMergeDiags.append(diagB)
						continue
			
				else:
					# endA[0] + currGap < startB[0]
				 	# stop the loop (remember that the diags are sorted!)
					# Impossible to merge next diags to diagA
					continueToSearch = False
					impossibleToMergeDiags.append(diagB)
			
			if len(fusionableDiags) > 0 :
				#assert all(x[4][0] >= endA[0] for x in fusionableDiags)
				#if dTa == '/':
				#	assert all( (x[0] == '/' or  x[0] == None ) for x in fusionableDiags)
				#	assert all(x[4][1]-endA[1] >= 0 for x in fusionableDiags), [(x[4][1],endA[1]) for x in fusionableDiags]
				#elif dTa == '\\':
				#	assert all( (x[0] == '\\' or x[0] == None ) for x in fusionableDiags)

				# Sort by priority depending on the choice of the distance :
				# either we want to fuse in priority along the diagonal line or relatively to the proximity on chromosome 1
				if distanceMetric in ['DPD', 'ED', 'CD', 'MD']:
					# score'Non'Diagonality = abs(deltaX - deltaY)
					scoreNonDiagonality= lambda x: abs(x[4][0]-diagA[5][0]) - abs(x[4][1]-diagA[5][1])
					fusionableDiags = sorted(fusionableDiags, key=scoreNonDiagonality, reverse=True)
				#sort by proximity to the diagonal line
				#elif distanceMetric == 'MD': # Doesn't change the result
				#	#sort by proximity on chromosome 1
				#	fusionableDiags = sorted(fusionableDiags, key=lambda x:x[4][0]-diagA[5][0])
				
				diagToFuse = fusionableDiags.pop(0) # Pop the first element of the list
				(dt,l1,l2,la,start,end) = diagToFuse
				# Lists are fused
				if dt!= None and dTa != None:
					if  dt == dTa:
						dt_res = dt
					else:
						# opposed diagTypes
						raise ErrorValue
				elif dt != None and dTa == None:
					dt_res = dt
				elif dTa != None and dt == None:
					dt_res = dTa
				else:
					# dt == None and dTa == None:
					# If the two diagonals are singletons
					dt_res = '/' if endA[1] < start[1] else '\\'
				la1.extend(l1)
				la2.extend(l2)
				laa.extend(la)
				endA = end
				diagA = (dt_res,la1,la2,laa,startA,endA)
				nbFusion += 1
				nbFusionCurrGap += 1
				#We still try to merge diags to the new diagA
				diagGen.putBack(diagA) # Needed to update diagA infos in the loop
			else:
				listOfFinishedDiags.append(diagA)
			# Diags are re-sorted owing to the A chromosome index
			for diag in sorted(fusionableDiags + impossibleToMergeDiags, key=lambda x:x[1][0]):
				diagGen.putBack(diag) # diagonals that were not fusionable are recorded to try to merge them after
			diagGen.rewind()
		print >> sys.stderr, "number of merges for currGap=%s :" % currGap, nbFusionCurrGap
	# Once all merges have been performed for a currGap it is necessary to repeat the merging process with all the diagonals for currGap+1 
	print >> sys.stderr, "Total number of merges", nbFusion
	print >> sys.stderr, "number of Diags after the merging process", len(sorted(list(diagGen) + listOfFinishedDiags, key=lambda x:x[1][0]))
	listOfSortedAndMergedDiagonals = sorted(list(diagGen) + listOfFinishedDiags, key=lambda x:x[1][0])
	return listOfSortedAndMergedDiagonals

def strandProduct(sa,sb):
		if sa != None and sb != None:
			return sa*sb
		else:
			return None

# Return the hm (or the mhp) of two compared chromosomes rewritten with family names instead of gene (or tb) names
# inputs :
# 	gcX : a chromosome of the genome X rewritten with the family names : [...,(f,s),...] 
# ouputs :
#	M : mh or mhp depending on whether gc1 and gc2 are chromosomes written in genes or in tbs
#		M is a dict of dict : M[i1][i2] = hpSign, this structure is lighter since M is a sparse matrix
#	locG2 : a dict {...,f:[(i2,s2)],...} with i2 the indices of the genes on gc2. For each f locG2 gives all its genes positions in the chromosome 2
# f : family, often the ancGene index
###################################################################################################################################
def homologyMatrix(gc1,gc2):
	locG2 = {}
	for (i2,(f,s2)) in enumerate(gc2):
        	if f != -1:
        		if f not in locG2:
        			locG2[f]=[]
        		locG2[f].append( (i2,s2) ) # LOCalisation on Gc 2
	#TODO optimization: Try to change the dict structure by a Matrix of len(gc1) rows and len(gc2) cols, maybe it will be faster
        M = {}
        for (i1,(f,s1)) in enumerate(gc1):
        	if f !=-1 and f in locG2: #TODO : remove by advance the f == -1 from gc1
        		for (i2,s2) in locG2[f]:
        			if i1 not in M:
        				M[i1]={}
        			M[i1][i2]= strandProduct(s1,s2)
	return (M, locG2) 

#
# At the beginning of a diagonal this function is called to set a diagonal type 
# inputs :
#	i1,i2 : the indices of the first hp of the current diagonal in M
#	M : the homology matrix (either mh or mhp)
#	consistentDiagonals : boolean (true if you want consistent diagonals otherwise false)
# output :
# 	diagType :
# 		either a bottom-left to top-right (slash : '/') diagonal 
# 		or top-left to bottom right (backslash : '\') diagonal
################################################################################################
def findDiagType(i1,i2,M,consistentSwDType):
	if M[i1][i2] != None and consistentSwDType:
		diagType = '/' if M[i1][i2] == +1 else '\\'
	else:
		if i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [1,None]) if consistentSwDType else True):
			diagType = '/'
		elif i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if consistentSwDType else True):
			diagType = '\\'
		else:
			diagType = None
	return diagType
#
# Extract sbs in a pairwise comparison of two chromosomes
############################################################
@utils.myTools.verbose
def extractSbsInPairCompChr(c1, c2, gc1, gc2, gapMax=0, distanceMetric = 'DPD', consistentSwDType=True, verbose = True):
	print >> sys.stderr, "(PPID = %s, PID = %s) start to extract diagonals on G1[%s]_vs_G2[%s]" % (os.getppid(), os.getpid(), c1, c2)
        listOfDiags = []
	(M,locG2) = homologyMatrix(gc1, gc2) 
        
	la=[]
        l1=[]
        l2=[]
        diagType = None
        strand = None
      	i1_old = -1
	
	#TODO scan M instead of gc1 : impossible since 'dict size cannot change during a loop over its items'
	# scan M from left to right
        for (i1,(f,_)) in enumerate(gc1):
		# f in locG2 means that i1 has at least one homolog on gc2, locG2 never changes, that is why we iterate over locG2
        	if f != -1 and f in locG2: 
			i1_old=i1
        		# scan M from bottom to top 
			for (i2,_) in locG2[f]: # the family name of the gene at index i2 is the same as f
				# When a diagonal is extracted the scanning of M must continue right after the first hp of the extracted diagonal
				i1=i1_old 
				# TODO write a recursive function easier understanding
        			# While hps can be added to a diagonal
        			while i1 in M and i2 in M[i1]: 
					# Here a diagonal is started or hps are added to an already started diagonal
					f = gc1[i1][0]
					if len(la) == 0:
						# first hp of a diagonal, the diagonalType is searched
						diagType = findDiagType(i1,i2,M,consistentSwDType)

					# A hp has a unique associated gene family (often a unique ancestral gene)
					# orientation of tbs on gc1 are used as a reference, if the orientation of a tb on gc1 is unknown, it is possible to infer the ancestral orientation by using the diagType and the orientation of the homologous tb on gc2
					if gc1[i1][1] != None:
						ancestralStrand = gc1[i1][1] 
					elif diagType == '/' and gc2[i2][1] != None:
						ancestralStrand = gc2[i2][1] 
					elif diagType == '\\' and gc2[i2][1] != None:
						ancestralStrand = -gc2[i2][1] 
					else:
						ancestralStrand = None

					la.append((f,ancestralStrand))
					l1.append((i1,gc1[i1][1]))
					l2.append((i2,gc2[i2][1]))
					del M[i1][i2]
					if len(M[i1].keys()) == 0:
						del M[i1] 
					#diagType == None if sameStrand == False and len(la) == 1, first element of a diagonal which we donnot take care of gene orientation
					#if (diagType == "/" or diagType == None) and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if consistentSwDType else True):
					if diagType == "/" and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if consistentSwDType else True):
						i1=i1+1
						i2=i2+1
						#assert i2-l2[-1][0] == +1
						#assert i1-l1[-1][0] == +1
						#assert i2 in M[i1]
					
					elif diagType=="\\" and i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if consistentSwDType else True):
						i1=i1+1
						i2=i2-1
						#assert i2-l2[-1][0] == -1
						#assert i1-l1[-1][0] == +1
						#assert i2 in M[i1]
					else: 
						# Since no more hps can be added to the current diagonal, the diagonal is recorded
						#assert len(la) > 0
						listOfDiags.append((diagType,l1,l2,la))
						l1=[]
						l2=[]
						la=[]
						diagType=None
						break # exit the while loop and iter the for loop
	# merging process, fuse diagonals 
	if len(listOfDiags) > 0 and gapMax >= 0 :
		listOfDiags = mergeSbs(listOfDiags, gapMax, distanceMetric=distanceMetric, verbose=verbose)
	listOfDiags = [((c1, diag[1]), (c2, diag[2]), diag[3]) for diag in listOfDiags] #FIXME diagType=diag[0] could be added as an information on the diagonal here
	return listOfDiags

#
# Wrapper used to print the % of the task progression
#######################################################
def wrapper_extractSbsPairCompChr(input, kwargs, output, NbOfTasks, listOfPercentage, lock):
	for args in iter(input.get, 'STOP'):
		result = extractSbsInPairCompChr(*args, **kwargs)
		output.put(result)
		avancement = 100-100*(input.qsize()) / NbOfTasks
		lock.acquire()
		if avancement in listOfPercentage:
			print >> sys.stderr, avancement, "% of the synteny block extraction"
			listOfPercentage.remove(avancement)
		lock.release()

#
# Multiprocess queue 
#####################
@utils.myTools.tictac
@utils.myTools.verbose
def extractSbsInPairCompGenomesMultiprocess(g1 , g2, gapMax=-1 ,distanceMetric='DPD', consistentSwDType=True, verbose=True):
	manager = Manager()
	lock = Lock() # Lock is used to give priority to one process in order that it is not disturbed by other process
	listOfPercentage = Manager().list() # Manager is used to manage objects which are shared among process
	for i in [10*j for j in range(11)]:
		listOfPercentage.append(i)

	NUMBER_OF_PROCESSES = multiprocessing.cpu_count() * 2
	TASKS = [(c1, c2, g1[c1], g2[c2]) for (c1,c2) in itertools.product([c1 for c1 in g1],[c2 for c2 in g2])]
	KWARGS = {'gapMax':gapMax, 'distanceMetric':distanceMetric, 'consistentSwDType':consistentSwDType, 'verbose':verbose}
	
	# Create queues
	task_queue = Queue()
	done_queue = Queue()
	
	# Submit tasks
	for task in TASKS :
		task_queue.put(task)

	# Start worker processes
	for i in range(NUMBER_OF_PROCESSES):
		Process(target=wrapper_extractSbsPairCompChr, args=(task_queue, KWARGS, done_queue, int(task_queue.qsize()), listOfPercentage, lock)).start()

	# Get and print results
	# 'Unordered results:'
	for i in range(len(TASKS)):
		for listOfDiags in  done_queue.get():
			yield listOfDiags

	# Tell child processes to stop
	for i in range(NUMBER_OF_PROCESSES):
		task_queue.put('STOP')
	
	return
	
# Rewrite genomes as a list of ancGeneIds
# ancGeneId can be considered as a family name
def rewriteWithAncGeneID(genome, ancGenes):
	newGenome = {}
	for c in genome.keys():
		#assert len(genome[c]) >=1
		# Genome.getPosition(name) returns a set of GenePosition that corresponds to the positions of the gene g.names in the ancGene file. 
		# Gene position is a namedtuple : GenePosition(chromosome=value, index=value)
		tmp = [(ancGenes.getPosition([gn]),s)  for (gn,s) in genome[c]]
		#assert all([len(g) <= 1 for (g,_) in tmp]) # DEBUG assertion, verify that there is only one ancGene ID for each gene
		#assert set(len(x[0]) for x in tmp).issubset(set([0,1]))
		#assert set(list(x[0])[0].chromosome for x in tmp if len(x[0]) > 0).issubset([None])
		newGenome[c] = [(g.pop().index if len(g) > 0 else -1, strand) for (g,strand) in tmp] # WARNING g.pop() remove and return a random element from set g
	return newGenome

# Depending on the filterType parameter:
#	- filterType = FilterType.None (genomes are not filtered)
#	- filterType = FilterType.InCommonAncestor (only genes herited from the ancestor are kept)
#	- filterType = FilterType.InBothSpecies (only 'anchor genes', ie genes present in both species, are kept)
# Returns (g1,g2,trans1,trans2)
# 	- g1 the genome 'g1' rewritten with ancGenes ID
#	- trans1 = { ..., newi:oldi, ...} with newi, the index of a gene in 'g1' and oldi the index of the same gene in the original genome 'g1'
def filter(g1_orig, g2_orig, filterType, minChromLength, keepOriginal=False):
	# Mark genes that are not in the intersection for future removal
	# Marking is done by switching (g,s) into (-1,s)
	def markGeneIntersection():
		def usedValues(genome): # create a set with all gene names
			val = set()
			for x in genome.itervalues():
				val.update(i for (i,_) in x)
			return val
		# Rewrite genomes using the intersection 
		# genes in intersection are kept other are changed into -1
		def rewrite(genome, inters): 
			for c in genome:
				genome[c] = [(i,s) if i in inters else (-1,s) for (i,s) in genome[c]]
		val1 = usedValues(g1)
		val2 = usedValues(g2)
		inters = val1.intersection(val2)
		rewrite(g1, inters)
		rewrite(g2, inters)
	# Remove too small chromosomes
	def filterSize(genome):
		flag = False
		for c in genome.keys():
			if len(genome[c]) < minChromLength:
				flag = True
				del genome[c]
		return flag
	# Remove genes that are not inherited of a gene in the ancestor
	# Warning : modify genome 
	def filterContent(genome):
		transNewToOld = collections.defaultdict(dict)
		for c in genome:
			tmp = [(i,x) for (i,x) in enumerate(genome[c]) if x[0] != -1] #i corresponds to indices of the non-filtered genome
			transNewToOld[c] = dict((newi,oldi) for (newi,(oldi,_)) in enumerate(tmp)) # here newi corresponds to indices of the filtered genome and oldi corresponds to indices of the non-filtered genome
			# oldi : index of gene before filtering
			# newi : index of gene after filtering
			genome[c] = [x for (_,x) in tmp]
		return transNewToOld
	
	if keepOriginal:
		g1 = copy.deepcopy(g1_orig)
		g2 = copy.deepcopy(g2_orig)
	else:
		g1=g1_orig
		g2=g2_orig
	g1filt2g1origin = collections.defaultdict(dict)
	g2filt2g2origin = collections.defaultdict(dict)

	# In all cases filtering on length
	# All genes are conserved
	if filterType == FilterType.None:
		filterSize(g1)
		filterSize(g2)
	# Only genes herited from an ancestral genes are conserved
	elif filterType == FilterType.InCommonAncestor:
		g1filt2g1origin = filterContent(g1)
		g2filt2g2origin = filterContent(g2)
		filterSize(g1)
		filterSize(g2)
	# Only conserve genes present in both extant genomes
	elif filterType == FilterType.InBothSpecies:
		while True:
			markGeneIntersection()
			g1filt2g1origin = filterContent(g1)
			g2filt2g2origin = filterContent(g2)
			useful = filterSize(g1) or filterSize(g2)
			if not useful: # if a chromosome has been removed because of the filtering on the length, the filtering is performed once more
				break
	else:
		# impossible case
		raise
	
	# In order to write genomes in files and verify sbs predictions
	# for c in g2:
	#	for (i,(g,s)) in enumerate(g2[c]):
	#		print >> sys.stderr, c, "\t", s, "\t", g2.lstGenes[c][trans2[c][i]].names[0] 

	return (g1, g2, g1filt2g1origin, g2filt2g2origin)

# Inputs :
#	genome_aID = [..., (aID,s), ...] with 'aID' the line number of the gene in the ancGene ans 's' the strand of the gene in the genome
# Outputs :
#	genome_tb and tb2g are two dics
#	for each chromosome c
#		genome_tb[c] = [..., (aID,s), ...] the rewritten genome into TBs. A TB is caracterized by its 'aID' and its strand, either 0,1 or None 
#		tb2g[c] (TB to genes) = [..., [i1,i2,i3], ...] contains as many lists as TBs. Each TB has a list [i1,i2,i3] which contains indices of the genes in 'genome_aID'
#		g2tb[c] (gene to TB) = [i_tb1, itb1, i_tb2, ...] for each gene of genome_aID the index of its tb in genome_tb
#####################################################################################################################################
def rewriteInTb(genome_aID):
	#TODO? it may be interesting to merge tbs that are closer or equal to a user defined trheshold 'tandemGap', introduce gappedTandemBlocks as in i-ADHoRe
	# in i-ADHoRe genes are remapped on the first paralog 
	
	# the rewritten genome
	genome_tb = {} # Need to keep dictionnaries because there is often a gap in chromosome notations
	tb2g = {}
	N_GTD_g = 0 # Number of Genic Tandem Duplication global
	for c in genome_aID:
		# temp values
		old_aID = None
		old_s = None
		genome_tb[c]=[]
		tb2g[c]=[]
		for i,(aID,s) in enumerate(genome_aID[c]):
			# if the current gene is not a paralog of the previous gene
			if aID != old_aID:
				# start a new TB
				genome_tb[c].append((aID,s))
				# add the current gene index to the new tb component of tb2g
				tb2g[c].append([i])
			else:
			# aID == old_aID:
				# if the current gene has a different orientation than the previous paralog, the orientation of the tb is set to unknown
				if s != old_s:
					genome_tb[c][-1] = (genome_tb[c][-1][0],None)
				# The index of the current gene is added to the last tb
				tb2g[c][-1].append(i)
				N_GTD_g += 1
			old_aID = aID
			old_s = s
	return (genome_tb, tb2g, N_GTD_g)

def convertGenomeIntoDicOfChromOfOrderedGenes(genome):
	newGenome = {} # Warning : it is important to use a dict since there are sometimes a jump in the numerotation of chromosomes in a genome
	for c in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:	
		assert len(genome.lstGenes[c]) >=1
		newGenome[c] = [(g.names[0],g.strand) for g in genome.lstGenes[c]]
	return newGenome

# compute the adviced maximal gap length (in genes) for the diagonal merger
# We compute the gap that gives the closest probability p-value(sb(nbHps,gap,N12,N1,N2)) to targetProba in an a MHP(N1,N2,N12). 
# Usually this MHP is an average of all the MHPs of the whole genome comparison. Even if considering this average MHP is not very rigorous, we only need reasonable rather than optimal values for recommendedGap.
@utils.myTools.verbose
def recommendedGap(nbHps, targetProba, N12, N1, N2, p_hpSign=None, maxGapThreshold=20, verbose=False):
	firstPrint=True
	tries = []
	for g in range(0,maxGapThreshold+1):
		L = (nbHps-1)*g+nbHps # lengths of the chromosomal windows
		pVal = utils.myProbas.pValue(nbHps,g,L,L,N12,N1,N2,p_hpSign,verbose=verbose)
		if firstPrint == True:
			print >> sys.stderr, "P-values of sbs of %s hps (nbHpsRecommendedGap) hps spaced by a gap 'g' in an average MHP(N1=%s, N2=%s, N12=%s) over all MHPs involved in the whole genome comparison.:" % (nbHps,N1,N2,N12)
			firstPrint=False
		print >> sys.stderr, "nbHps=%s, g=%g, pVal=%s" % (nbHps, g,pVal)
		if pVal != None:
			tries.append(pVal)
		elif pVal == None:
			tries.append(-float('inf'))
		else:
			raise ValueError
	print >> sys.stderr, "The first p-value over targetProba=%s (targetProbaRecommendedGap) sets the recommended gap value" % targetProba
	for (g,pVal) in  enumerate(tries):
		if pVal > targetProba:
			return g
	#min_index, min_value = min(enumerate(tries), key=lambda x:x[1])
	#g = min_index
	#return g

# Number of non-empty value in the mh or the mhp
@utils.myTools.tictac
def numberOfHomologies(g1,g2):
	nbHomologies = {}
	for c1 in g1:
		for c2 in g2:
			comp = (c1,c2)
			gc1 = g1[c1]
			gc2 = g2[c2]
			(Ms,_) = homologyMatrix(gc1, gc2)
			nbHomologies[comp] = sum([len(Ms[i1]) for i1 in Ms])
	nbHomologies_g = sum([a for a in nbHomologies.values()])
	#assert nbHomologies_g >= min(sum([len(g1[c]) for c in g1]), sum([len(g2[c]) for c in g2])),"%s,%s" %  (sum([len(g1[c]) for c in g1]), sum([len(g2[c]) for c in g2])) 
	#Not needed since the two lineages may have undergone some differential gene losses
	return nbHomologies, nbHomologies_g

#returns the maximum gap between two hps with the Chebyschev Distance metric (CD)
def max_g_of(diag):
	if len(diag) == 3:
		((c1,l1_tb),(c2,l2_tb),la_tb) = diag
	elif len(diag) == 4:
		((c1,l1_tb),(c2,l2_tb),la_tb,p) = diag
	l1s = [(l1_tb[i1][0],l1_tb[i1+1][0]) for i1,_ in enumerate(l1_tb[:-1])]
	max_g1 = max(l1s,key=lambda x:abs(x[1]-x[0]))
	max_g1 = abs(max_g1[0]-max_g1[1]) - 1 #do not forget the -1
	l2s = [(l2_tb[i2][0],l2_tb[i2+1][0]) for i2,_ in enumerate(l2_tb[:-1])]
	max_g2 = max(l2s,key=lambda x:abs(x[1]-x[0]))
	max_g2 = abs(max_g2[0]-max_g2[1]) - 1 #do not forget the -1
	max_g = max(max_g1,max_g2)
	return max_g

# Statistical validation of synteny blocks
# because of collections.Counter
@utils.myTools.minimalPythonVersion((2,7))
@utils.myTools.verbose
def filterStatisticalValidation(listOfDiags, g1_tb, g2_tb, N12s, p_hpSign, pThreshold=0.001, NbOfHomologiesThreshold=50, validateImpossToCalc_mThreshold=3, verbose=True):
	def  calculateCharacteristics(diag):
		if len(diag) == 3:
			((c1,l1_tb),(c2,l2_tb),la_tb) = diag
		elif len(diag) == 4:
			((c1,l1_tb),(c2,l2_tb),la_tb,p) = diag
		comp = (c1,c2)
		m = len(la_tb) # number of homologies in the sb
		#(x0,y0) = (l1_tb[0][0],l2_tb[0][0])
		#(x1,y1) = (l1_tb[1][0],l2_tb[1][0])
		#diagType = '/' if (x1 - x0)*(y1 - y0) else '\\'
		lw1=[sys.maxint,-sys.maxint]
		lw2=[sys.maxint,-sys.maxint]
		for tb1 in l1_tb:
			lw1[0] = tb1[0] if  tb1[0] < lw1[0] else lw1[0]
			lw1[1] = tb1[0] if  tb1[0] > lw1[1] else lw1[1]
		for tb2 in l2_tb:
			lw2[0] = tb2[0] if  tb2[0] < lw2[0] else lw2[0]
			lw2[1] = tb2[0] if  tb2[0] > lw2[1] else lw2[1]
		lw1 = lw1[1]-lw1[0]+1 # +1 in order to have the nb of genes in the W1
		lw2 = lw2[1]-lw2[0]+1 # +1 --------------------------------------- W2
		assert lw1 >= 0, lw1
		assert lw2 >= 0, lw2
		l1_min=min(l1_tb,key=lambda x:x[0])[0]
		l1_max=max(l1_tb,key=lambda x:x[0])[0]
		l2_min=min(l2_tb,key=lambda x:x[0])[0]
		l2_max=max(l2_tb,key=lambda x:x[0])[0]
		if m > 1:
			assert l1_tb[1] >=  l1_tb[0] # Because of the scanning process
			max_g = max_g_of(diag)
		else:
			max_g = None
		return (m,max_g,lw1,lw2,l1_min,l1_max,l2_min,l2_max)

	listOfRejectedDiags = []
	listOfImpossibleToCalcProbaDiags = []
	listOfStatValDiags = []
	for i,diag in enumerate(listOfDiags):
		((c1,l1_tb),(c2,l2_tb),la_tb) = diag
		(m,max_g,lw1,lw2,l1_min,l1_max,l2_min,l2_max) = calculateCharacteristics(diag)
		comp = (c1,c2)
		na = len(g1_tb[c1]) # Nb of Tbs on C1
		nb = len(g2_tb[c2]) # Nb of Tbs on C2
		nab = N12s[comp] # Nb of homologies in the mhp	
		if m <= 1:
			# all diagonals of length 1 are rejected 
			listOfRejectedDiags.append((diag[0],diag[1],diag[2],None))
			continue 
		elif m > nab:
			# there are not m hps in the MHP
			listOfImpossibleToCalcProbaDiags.append((diag[0],diag[1],diag[2],None))
			continue
		elif m > min([lw1, lw2]):
			# there are too many dispersed paralogies in the window
			listOfImpossibleToCalcProbaDiags.append((diag[0],diag[1],diag[2],None))
			continue

		if m > NbOfHomologiesThreshold: # This is to avoid too time consuming computations 
			p=0
		else:
			comp = (c1,c2)
			p = utils.myProbas.pValue(m, max_g, lw1, lw2, nab, na, nb, p_hpSign[comp], verbose=True)
		if p == None:
			listOfImpossibleToCalcProbaDiags.append((diag[0],diag[1],diag[2],None))
		elif p < pThreshold:
			listOfStatValDiags.append((diag[0],diag[1],diag[2],p))
		else:
			listOfRejectedDiags.append((diag[0],diag[1],diag[2],p))
			assert len(diag[2]) == len(diag[0][1])
			assert len(diag[2]) == len(diag[1][1])

	# automatically validate sbs that contain more than validateImpossToCalc_mThreshold hps
	listOfImpossibleToCalcProbaDiagsStatVal=[]
	listOfImpossibleToCalcProbaDiagsRejected=[]
	for diag in listOfImpossibleToCalcProbaDiags:
		m = calculateCharacteristics(diag)[0]
		if m >= validateImpossToCalc_mThreshold: 
			# if the diagonal contains more than validateImpossToCalc_mThreshold  hps, it is validated in order to avoid to reject long and perfect diagonals because of only one dispersed tandem duplication
			listOfImpossibleToCalcProbaDiagsStatVal.append(diag)
		else:
			listOfImpossibleToCalcProbaDiagsRejected.append(diag)
	
	# update the lists
	listOfStatValDiags = listOfStatValDiags + listOfImpossibleToCalcProbaDiagsStatVal
	listOfRejectedDiags = listOfRejectedDiags + listOfImpossibleToCalcProbaDiagsRejected

	assert len(listOfStatValDiags) + len(listOfRejectedDiags) == len(listOfDiags)
	assert len(listOfImpossibleToCalcProbaDiags) == len(listOfImpossibleToCalcProbaDiagsStatVal) + len(listOfImpossibleToCalcProbaDiagsRejected)

	firstPrint = True
	for diag in listOfRejectedDiags:
		((c1,l1_tb),(c2,l2_tb),la_tb,p) = diag
		(m,max_g,lw1,lw2,l1_min,l1_max,l2_min,l2_max) = calculateCharacteristics(diag)
		comp = (c1,c2)
		na = len(g1_tb[c1]) # Nb of Tbs on C1
		nb = len(g2_tb[c2]) # Nb of Tbs on C2
		nab = N12s[comp] # Nb of homologies in the mhp	
		if m>1:
			if firstPrint:
				print >> sys.stderr, "Diagonals of more than 1 hp that have been rejected during the statistical test :"
				firstPrint = False
			# TODO comment se fait-il que certaines diagonales soient en double ici ?
			print >> sys.stderr, "(c1=%s:%s-%s,c2=%s:%s-%s) \t (m=%s, max_g=%s, lw1=%s lw2=%s, nab=%s, na=%s, nb=%s)" % (c1,l1_min,l1_max,c2,l2_min,l2_max,m,max_g,lw1,lw2,nab,na,nb), "\t p=", p
	
	firstPrint = True
	for diag in listOfStatValDiags:
		((c1,l1_tb),(c2,l2_tb),la_tb,p) = diag
		(m,max_g,lw1,lw2,l1_min,l1_max,l2_min,l2_max) = calculateCharacteristics(diag)
		if m == 2:
			comp = (c1,c2)
			na = len(g1_tb[c1]) # Nb of Tbs on C1
			nb = len(g2_tb[c2]) # Nb of Tbs on C2
			nab = N12s[comp] # Nb of homologies in the mhp		
			if firstPrint:
				print >> sys.stderr, "Diagonals containing 2 hps that have passed the statistical test :"
				firstPrint = False
			print >> sys.stderr, "(c1=%s:%s-%s,c2=%s:%s-%s) \t (m=%s, max_g=%s, lw1=%s lw2=%s, nab=%s, na=%s, nb=%s)" % (c1,l1_min,l1_max,c2,l2_min,l2_max,m,max_g,lw1,lw2,nab,na,nb), "\t p=", p

	def statsDiagsM(listOfDiagsX,m):
		if listOfDiagsX == []:
			return []
		if len(listOfDiagsX[0]) == 4:
			listOfDiagsXM = [((c1,l1_tb),(c2,l2_tb),la_tb) for ((c1,l1_tb),(c2,l2_tb),la_tb,p) in listOfDiagsX if len(la_tb)==m]
		elif len(listOfDiagsX[0]) == 3:
			listOfDiagsXM = [((c1,l1_tb),(c2,l2_tb),la_tb) for ((c1,l1_tb),(c2,l2_tb),la_tb) in listOfDiagsX if len(la_tb)==m]
		listOfDiagsXMGaps = [max_g_of(diag) for diag in listOfDiagsXM]
		diagsXMGaps = collections.Counter(listOfDiagsXMGaps)
		diagsXMGaps = ["%s:%s" % (length,nb) for (length,nb) in sorted(diagsXMGaps.items())]
		return diagsXMGaps
	print >> sys.stderr, "Over all sbs of 2 hps before stat. val. the distribution of gap maximum is : {%s}" % " ".join(statsDiagsM(listOfDiags,2))
	print >> sys.stderr, "Over all rejected sbs of 2 hps the distribution of gap maximum is : {%s}" % " ".join(statsDiagsM(listOfRejectedDiags,2))
	print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps the distribution of gap maximum is : {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(listOfImpossibleToCalcProbaDiags,2))
	print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps which are finally validated, the distribution of gap maximum is : {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(listOfImpossibleToCalcProbaDiagsStatVal,2))
	print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps which are finally rejected, the distribution of gap maximum is : {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(listOfImpossibleToCalcProbaDiagsRejected,2))

	print >> sys.stderr, "Over all stat. val. sbs of 2 hps the distribution of gap maximum is : {%s}" % " ".join(statsDiagsM(listOfStatValDiags,2))
	
	def statsDiagLengths(listOfDiagsX):
		if listOfDiagsX == []:
			return []
		if len(listOfDiagsX[0]) == 4:
			listOfDiagsX = [((c1,l1_tb),(c2,l2_tb),la_tb) for ((c1,l1_tb),(c2,l2_tb),la_tb,p) in listOfDiagsX]
		elif len(listOfDiagsX[0]) == 3:
			listOfDiagsX = [((c1,l1_tb),(c2,l2_tb),la_tb) for ((c1,l1_tb),(c2,l2_tb),la_tb) in listOfDiagsX]
		diagXLengths = collections.Counter([len(la_tb) for (_,_,la_tb) in listOfDiagsX])
		diagXLengths = ["%s:%s" % (length,nb) for (length,nb) in sorted(diagXLengths.items())]
		return diagXLengths
	print >> sys.stderr, "Over all diagonal before the stat. val., distribution of the diag lengths : {%s}" % " ".join(statsDiagLengths(listOfDiags))
	print >> sys.stderr, "Over all rejected diagonals, distribution of all diag lengths : {%s}" %  " ".join(statsDiagLengths(listOfRejectedDiags))
	print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies), distribution of diag lengths : {%s}" %  " ".join(statsDiagLengths(listOfImpossibleToCalcProbaDiags))
	print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies) which are finally validated, distribution of diag lengths : {%s}" %  " ".join(statsDiagLengths(listOfImpossibleToCalcProbaDiagsStatVal))
	print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies) which are finally rejected, distribution of diag lengths : {%s}" %  " ".join(statsDiagLengths(listOfImpossibleToCalcProbaDiagsRejected))
	print >> sys.stderr, "Over all diagonals that passed the stat. val., distribution of diag lengths : {%s}" %  " ".join(statsDiagLengths(listOfStatValDiags))

	return listOfStatValDiags

# because of collections.Counter
@utils.myTools.minimalPythonVersion((2,7))
def numberOfDuplicates(g_aID_filt):
	# create a long chromosome with all the chromosome genes
	newG=[]
	for chr in g_aID_filt.values():
		chr =  [g for g,_ in chr]
		newG = newG + chr
	N_GD_g = sum([duplicates-1 for duplicates in collections.Counter(newG).values()])
	return (N_GD_g)

#
# Complete procedure to compute synteny blocks (sbs) between 2 genomes, using homology relationships contained in ancGenes
# Inputs:
#	g1, g2 : utils.myGenomes.Genome (or dict : {...c:[..., (g,s), ...]}) of the two compared species
#	ancGenes : utils.myGenomes.Genome define the homology relationships (usually the ancGenes of the LCA of the two compared species)
#	consistentSwDType :	True  => gene transcription orientations must be consistent with diagonal types ('/' or '\')
#				False => Do not take care of gene orientation
#	gapMax : distance (in tandemblocks) allowed between two diagonals to be merged together. This allows to go over wrong  annotations. But it can also introduce some Fp
# 	minChromLength is the minimal length of the chromosome to be considered
#	distanceMetric : the distance metric (either MD,DPD,ED or CD)
#	pThreshold : the probability threshold under which the sb is considered significant
#	FilterType 
#		InCommonAncestor : all genes herited from a gene of the LCA are kept during the filtering of extant genomes
#		InBothGenomes : only 'anchor genes' (genes present in both genomes) are kept
#		None : genomes are not filtered
# Outputs:
#	synteny blocks are yielded as (strand,(c1,l1),(c2,l2),la)
#	c1 : chromosomes of genome 'g1' where there is a synteny block corresponding to the diagonal
#	l1 : the synteny block on genome 'g1' 
#	l1 = [...,(i,s),...] with 'i' the index of the gene and 's' its strand in the genome 1 without species specific genes
#
#########################################################################################################################
@utils.myTools.tictac
@utils.myTools.verbose
def extractSbInPairCompGenomes(g1, g2, ancGenes, gapMax=None, distanceMetric='DPD', pThreshold=0.001, filterType=FilterType.None, consistentSwDType=True, minChromLength=0, nbHpsRecommendedGap=2, targetProbaRecommendedGap=0.01, validateImpossToCalc_mThreshold=3, multiProcess=True, verbose=True):

	if isinstance(g1,utils.myGenomes.Genome) and isinstance(g2,utils.myGenomes.Genome):
		g1 = convertGenomeIntoDicOfChromOfOrderedGenes(g1)
		g2 = convertGenomeIntoDicOfChromOfOrderedGenes(g2)
	elif isinstance(g1,dict) and isinstance(g2,dict):
		pass
	else:
		raise TypeError('g1 and/or g2 must be either utils.myGenomes.Genome or dict')
	
	#step 1 : remove genes specific to each lineage
	################################################
	# rewrite genomes by family names (ie ancGene names)
	g1_aID = rewriteWithAncGeneID(g1, ancGenes) 
	g2_aID = rewriteWithAncGeneID(g2, ancGenes)
	# genes that are not in ancGene have a aID=-1
	print >> sys.stderr, "genome 1 initially contains %s genes" %  sum([len(g1[c1]) for c1 in g1])
	print >> sys.stderr, "genome 2 initially contains %s genes" %  sum([len(g2[c2]) for c2 in g2])
	(g1_aID_filt, g2_aID_filt, filt2origin1, filt2origin2) = filter(g1_aID, g2_aID, filterType, minChromLength) # Must be applied on the two genomes, because of the mode inBothGenomes (InCommonAncestor => not only anchor genes are kept but all genes herited from a gene of the LCA)
	print >> sys.stderr, "genome 1 after filterType=%s contains %s genes" % (filterType, sum([len(g1_aID_filt[c1]) for c1 in g1_aID_filt]))
	print >> sys.stderr, "genome 2 apres filterType=%s contains %s genes" % (filterType, sum([len(g2_aID_filt[c2]) for c2 in g2_aID_filt]))
	N_GD_1_g  = numberOfDuplicates(g1_aID_filt)
	N_GD_2_g  = numberOfDuplicates(g2_aID_filt)
	print >> sys.stderr, "genome 1 contains %s gene duplicates (initial gene excluded)" % N_GD_1_g
	print >> sys.stderr, "genome 2 contains %s gene duplicates (initial gene excluded)" % N_GD_2_g
	(g1_tb, tb2g1, N_GTD_1_g) = rewriteInTb(g1_aID_filt)
	(g2_tb, tb2g2, N_GTD_2_g) = rewriteInTb(g2_aID_filt)
	print >> sys.stderr, "genome 1 rewritten in tbs, contains %s tbs" % sum([len(g1_tb[c1]) for c1 in g1_tb])
	print >> sys.stderr, "genome 2 rewritten in tbs, contains %s tbs" % sum([len(g2_tb[c2]) for c2 in g2_tb])
	N12s, N12_g = numberOfHomologies(g1_tb,g2_tb)
	print >> sys.stderr, "pairwise comparison of genome 1 and genome 2 yields %s hps" % N12_g
	print >> sys.stderr, "genome 1 contains %s tandem duplicated genes (initial gene excluded)" % N_GTD_1_g
	print >> sys.stderr, "genome 2 contains %s tandem duplicated genes (initial gene excluded)" % N_GTD_2_g
	N_DD_1_g  = numberOfDuplicates(g1_tb)
	N_DD_2_g  = numberOfDuplicates(g2_tb)
	print >> sys.stderr, "genome 1 contains %s dispersed duplicated tbs (initial tb excluded)" % N_DD_1_g
	print >> sys.stderr, "genome 2 contains %s dispersed duplicated tbs (initial tb excluded)" % N_DD_2_g
	assert N_DD_1_g + N_GTD_1_g == N_GD_1_g
	assert N_DD_2_g + N_GTD_2_g == N_GD_2_g

	# compute the recommended gapMax parameter
	#########################################
	(p_hpSign,p_hpSign_g,(sTBG1, sTBG1_g),(sTBG2, sTBG2_g)) = utils.myProbas.statsHpSign(g1_tb,g2_tb)
	print >> sys.stderr, "genome 1 tb orientation proba = {+1:%s,-1:%s,None:%s} (stats are also calculated for each chromosome)" % (sTBG1_g[+1], sTBG1_g[-1], sTBG1_g[None])
	print >> sys.stderr, "genome 2 tb orientation proba = {+1=%s,-1:%s,None:%s} (stats are also calculated for each chromosome)" % (sTBG2_g[+1], sTBG2_g[-1], sTBG2_g[None])
	print >> sys.stderr, "hp sign proba in the 'global' mhp = {+1:%s,-1:%s,None:%s) (probabilities are also calculated for pairwise mhp)" % (p_hpSign_g[+1], p_hpSign_g[-1], p_hpSign_g[None])
	N1_g=sum([len(g1_tb[c1]) for c1 in g1_tb])
	N2_g=sum([len(g2_tb[c2]) for c2 in g2_tb])
	#Build an average MHP
	#N50 is better than the mean of the chromosome lengths, since it is less sensitive to numerous and very small chromosome lengths
	#N1_N50 = utils.myMaths.myStats.getValueNX(sorted([len(g1_tb[c1]) for c1 in g1_tb]),50) # calculate the N50
	#N2_N50 = utils.myMaths.myStats.getValueNX(sorted([len(g2_tb[c2]) for c2 in g2_tb]),50) # calculate the N50
	#weightedAverage is even better than N50 since it returns a more stable value, not a length of a chromosome of the karyotype, it better reflects the global distribution
	def weightedAverage(listOfChrLengths):
		sumLengths = sum(listOfChrLengths)
		return int(sum([(float(chrLength)/sumLengths)*chrLength for chrLength in listOfChrLengths]))

	#Waring: if the genome is badly assembled and contain a lot of small contigs, this averaging is not relevant and the recommended gapMax won't be relevant neither
	if len(g1_tb) > 50 or len(g2_tb) > 50:
		print >> sys.stderr, "Warning: one of the two genome seems badly assembled, this may mislead the recommended gapMax calculation"
	N1_weightedAverage = weightedAverage([len(g1_tb[c1]) for c1 in g1_tb])
	N2_weightedAverage = weightedAverage([len(g2_tb[c2]) for c2 in g2_tb])
	density = float(N12_g)/(N1_g*N2_g)
	N12_weightedAverage = int(density*N1_weightedAverage*N2_weightedAverage) # conservation of the density
	gap = recommendedGap(nbHpsRecommendedGap, targetProbaRecommendedGap, N12_weightedAverage, N1_weightedAverage, N2_weightedAverage, p_hpSign=p_hpSign_g, verbose=verbose)
	print >> sys.stderr, "recommended gapMax = %s tbs" % gap
	if gapMax == None:
		gapMax = gap
	print >> sys.stderr, "used gapMax = %s tbs" % gapMax
	
	# step 2 and 3 : build the MHP and extract putative sbs as diagonals
	#################################################################################
	# extract sbs in the tb base
	listOfSbs = []
	if multiProcess and (len(g1_tb.keys()) > 1 or len(g2_tb.keys()) > 1):
		# if the multiprocess option is True and if there is more than one pairwise comparison of chromosomes
		for tmpListOfSbs in extractSbsInPairCompGenomesMultiprocess(g1_tb, g2_tb, gapMax=gapMax, distanceMetric=distanceMetric, consistentSwDType=consistentSwDType, verbose=False):
			listOfSbs.append(tmpListOfSbs)
	else:
		for (chr1,chr2) in itertools.product(g1_tb.keys(), g2_tb.keys()):
			tmpListOfSbs = extractSbsInPairCompChr(chr1, chr2, g1_tb[chr1], g2_tb[chr2], gapMax=gapMax, distanceMetric=distanceMetric, consistentSwDType=consistentSwDType, verbose=False)
			listOfSbs.extend(tmpListOfSbs)

	# setp 4 : statistical validation of putative sbs
	##################################################
	#DEBUG
	#print >> sys.stderr, "sTBG1['Y']=%s" % sTBG1['Y']
	#print >> sys.stderr, "sTBG2['Y']=%s" % sTBG2['Y']
	#print >> sys.stderr, "p_hpSign[('Y','Y')]=%s" % p_hpSign[('Y','Y')]
	sbsGenSV = filterStatisticalValidation(listOfSbs, g1_tb, g2_tb, N12s, p_hpSign, pThreshold = pThreshold, NbOfHomologiesThreshold=50, validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold, verbose=verbose) # SV : statistical validation

	# format output
	for (i,diag) in enumerate(sbsGenSV):
		# translate from the tb base to gene base
		((c1,l1_tb),(c2,l2_tb),la_tb,pVal) = diag
		tb2gc1 = tb2g1[c1]
		tb2gc2 = tb2g2[c2]
		gc1_aID_filt = g1_aID_filt[c1]
		gc2_aID_filt = g2_aID_filt[c2]
		filt2originc1 = filt2origin1[c1]
		filt2originc2 = filt2origin2[c2]
		la=[]
		l1=[]
		l2=[]
		# convert diags from tb to genes before yielding them
		# each hp corresponds to an ancestral gene
		for (indx_HP,aGene) in enumerate(la_tb):
			(indx_tb_g1,_) = l1_tb[indx_HP]
			l1.extend([(indx,gc1_aID_filt[indx][1]) for indx in tb2gc1[indx_tb_g1]])
			(indx_tb_g2,_) = l2_tb[indx_HP]
			l2.extend([(indx,gc2_aID_filt[indx][1]) for indx in tb2gc2[indx_tb_g2]])
			la.append((ancGenes.lstGenes[None][aGene[0]].names[0], aGene[1], len(tb2gc1[indx_tb_g1]), len(tb2gc2[indx_tb_g2]))) # add informations on the size of the HP.
		if filterType != FilterType.None:
			sbsGenSV[i] = ((c1,[(filt2originc1[i1],s1) for (i1,s1) in l1]), (c2,[(filt2originc2[i2],s2) for (i2,s2) in l2]), la, pVal)
		else:
			sbsGenSV[i] = ((c1,l1), (c2,l2), la, pVal)
	return sbsGenSV

#################################
# Post processing of diags 
#################################

# Build a genome with sbs as contigs and only conserve one ancestral gene by hp
def buildGenomeFromSbs(listOfSbs, sbLengthThreshold):
	cptSb=0
	for sb in listOfSbs:
		((c1,l1),(c2,l2),la)
		old_ancGene=None
		for ihp,ancGene in enumerate(la):
			#assert ancGene != old_ancGene, "%s,%s" % (ancGene, old_ancGene) # Happens when there are paralogous hps in the same diagonal
			if ancGene == old_ancGene:
				continue # This may create diags of size 1 that will need to be removed
			ancSbGenome[cptSb].append((ancGene,ss[ihp])) # ancGene synonymous of hp-name
			old_ancGene = ancGene
		
		if len(ancSbGenome[cptSb]) <= 1: # In every cases sbs of length 1 are removed
			del ancSbGenome[cptSb]
			continue
		elif sbLengthThreshold != None and len(ancSbGenome[cptSb]) <= sbLengthThreshold:
			del ancSbGenome[cptSb]
			continue
		cptSb = cptSb+1
	return ancSbGenome

@utils.myTools.verbose
def computeAncestralCoverageBySbs(g1_tb, g2_tb, ancSbGenome, verbose = True):
	ancGenesInSbs=set([])
	for ancSb in ancSbGenome.values():
		for (ancGene,_) in ancSb:
			ancGenesInSbs.add(ancGene)
	print >> sys.stderr, "Nb of ancGenes in synteny blocks (each ancGene can appear at most once)", len(ancGenesInSbs)
	ancGenesInG1Tb = set([])
	ancGenesInG2Tb = set([])
	for tb1 in [tb1 for c1 in g1_tb for tb1 in g1_tb[c1]]:
		(ancGene,s) = tb1
		ancGenesInG1Tb.add(ancGene)
	for tb2 in [tb2 for c2 in g2_tb for tb2 in g2_tb[c2]]:
		(ancGene,s) = tb2
		ancGenesInG2Tb.add(ancGene)
	#{ancGenes that are present in G1 rewritten in tbs and in G2 rewritten in tbs}
	ancGenesInG1TbAndInG2Tb = ancGenesInG1Tb & ancGenesInG2Tb # intersection
	print >> sys.stderr, "coverage LCA_S1-S2 = cardinal({ancGenes in Sb}) / cardinal({ancGenes that are present in G1 rewritten in tbs and in G2 rewritten in tbs})\ncoverage LCA_S1-S2 = ", float(len(ancGenesInSbs)) / len(ancGenesInG1TbAndInG2Tb)
	coverage = float(len(ancGenesInSbs)) / len(ancGenesInG1TbAndInG2Tb)
	return coverage
