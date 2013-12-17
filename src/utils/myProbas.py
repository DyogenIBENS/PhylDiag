# -*- coding: utf-8 -*-
# PhylDiag v1.0
# needs python 2.7 at least
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software; you may copy, modify and/or distribute this work
# under the terms of the GNU General Public License, version 2 or later

######################
# probability package
######################

##################################
# Probability for synteny blocks
##################################
# This code is based on the publication ``PhylDiag : identifying complex cases of conserved synteny that include tandem duplications''
# authors :
# Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
# Matthieu MUFFATO (EBI, Cambridge, United Kingdom, muffato@ebi.ac.uk)
# Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)

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

from utils.myMaths import combinations
import collections
import math

# p_tbO_g={+1:%+1, -1:%-1, None:%None}, percentages of each tb orientation in genome_tb
# p_tbO is a dict containing the percentages for each chromosome
def statsTbOrientation(genome_tb):
	p_tbO = {}
	nbTb_plus_g=0
	nbTb_minus_g=0
	nbTb_None_g=0
	for c in genome_tb:
		nbTb_plus=0
		nbTb_minus=0
		nbTb_None=0
		for tb in genome_tb[c]:
			if tb[1] == +1:
				nbTb_plus += 1
			elif tb[1] == -1:
				nbTb_minus += 1
			elif tb[1] == None:
				nbTb_None += 1
			else:
				raise ValueError()
		assert len(genome_tb[c]) == nbTb_plus + nbTb_minus + nbTb_None
		p_tbO[c] = [float(v)/len(genome_tb[c]) for v in (nbTb_plus,nbTb_minus,nbTb_None)]
		p_tbO[c] = dict(  ((+1,p_tbO[c][0]),(-1,p_tbO[c][1]),(None,p_tbO[c][2])) ) 
		nbTb_plus_g += nbTb_plus
		nbTb_minus_g += nbTb_minus
		nbTb_None_g += nbTb_None
	nbTbs_g = sum([len(genome_tb[c]) for c in genome_tb])
	assert nbTbs_g == nbTb_plus_g + nbTb_minus_g + nbTb_None_g
	p_tbO_g = [float(v)/nbTbs_g for v in (nbTb_plus_g,nbTb_minus_g,nbTb_None_g)]
	p_tbO_g = dict(  ((+1,p_tbO_g[0]),(-1,p_tbO_g[1]),(None,p_tbO_g[2])) )
	assert all(abs(sum(p_tbO[c].values()) - 1) < 0.001 for c in p_tbO)
	assert abs(sum(p_tbO_g.values())-1) < 0.001
	return (p_tbO, p_tbO_g)

# returns p_hpSign_g {+1:%+1, -1:%-1, None:%None}, percentages of each hp sign in the 'global' mhp of the pairwise comparison of the two genomes
# p_hpSign is a dict containing the stats for each pairwise comparisons of chromosomes
def statsHpSign(g1_tb,g2_tb):
	# compute statistics on tandem blocks orientations
	sTbG1, sTbG1_g = statsTbOrientation(g1_tb)
	sTbG2, sTbG2_g = statsTbOrientation(g2_tb)
	# Probability of hps signs
	p_hpSign = collections.defaultdict(dict)
	for c1 in g1_tb:
		for c2 in g2_tb:
			comp = (c1,c2)
			# Here we can use sets since a hp sign +1 in a comparison between ca on X vs cb on Y is still a hp sign +1 in a comparison between cb on X vs ca on Y 
			p_hpSign[comp][+1] = sTbG1[c1][+1]*sTbG2[c2][+1] + sTbG1[c1][-1]*sTbG2[c2][-1]
			p_hpSign[comp][-1] = sTbG1[c1][+1]*sTbG2[c2][-1] + sTbG1[c1][-1]*sTbG2[c2][+1]
			p_hpSign[comp][None] = 	sTbG1[c1][+1]*sTbG2[c2][None] +  sTbG1[c1][None]*sTbG2[c2][+1] + \
					sTbG1[c1][-1]*sTbG2[c2][None] +  sTbG1[c1][None]*sTbG2[c2][-1] + \
					sTbG1[c1][None]*sTbG2[c2][None]
	p_hpSign_g={}
	for hpSign in [+1,-1,None]:
		p_hpSign_g[hpSign] = float(sum([p_hpSign[comp][hpSign] for comp in p_hpSign]))/len(p_hpSign)
	assert all(abs(sum(p_hpSign[comp].values()) - 1) < 0.001 for comp in p_hpSign)
	return (p_hpSign,p_hpSign_g,(sTbG1, sTbG1_g),(sTbG2, sTbG2_g))

# probability that there are exactly k hps in the window Wab of width la and height lb in a mhp of width na height nb and containing nab hps
# Warning : an assumption is made that in the mhp there are only orthologies and no paralogies
def p_d(k,la,lb,nab,na,nb):
	if k > nab:
		return None # In case of paralogies
	limSum = min(la-k,nab-k)
	sum_ = 0
	for i in range(0,limSum+1):
		try:
			foo = combinations(nab-k,i)*combinations(na-nab,la-(k+i))*combinations(nb-(k+i),lb-k)
		except:
			print >> sys.stderr, "Warning : not able to compute p_d(%s,%s,%s,%s,%s,%s)" % (k,la,lb,nab,na,nb)
			foo = 0
		sum_ += foo

	num = combinations(nab,k)*sum_ 
	den = combinations(na,la)*combinations(nb,lb)
	if den == float('inf') and num != float('inf'):
		return 0.0
	elif den == float('inf') and num == float('inf'):
		return None
	elif den == 0.0 and num != float('inf'):
		return float('inf')
	elif den == 0.0 and num == float('inf'):
		return None
	else:
		return float(num)/den

d_0 = lambda k,g,l : sum([((-1)**i)*combinations(k-1,i)*combinations(l-i*(g+1),k) for i in range(0,int(math.floor(float(l-k)/float(g+1))+1))])

# probability that k marked tbs form a maxgap(g)-cluster in a sequence of l tbs
def p_g_1D(k,g,l):
	if k > l: # if there are paralogies
		return None
	w = (k-1)*g+k
	num = max(0,l-w+1)*math.pow((g+1),k-1)+d_0(k,g,min(l,w-1))
	den = combinations(l,k)
	return float(num)/den

# probability that k marked hps form a maxgap(g)-cluster in both sequences of lengths la and lb
def p_g_2D(k,g,la,lb):
	try:
		return p_g_1D(k,g,la) * p_g_1D(k,g,lb)
	except TypeError:
		return None

# probability that k hps that are close enough form a slash diagonal
def p_slash(k,p_hpSign):
	return 1.0/math.factorial(k)*math.pow(p_hpSign[+1]+p_hpSign[None],k)
# probability that k hps that are close enough form a backslash diagonal
def p_backslash(k,p_hpSign):
	return 1.0/math.factorial(k)*math.pow(p_hpSign[-1]+p_hpSign[None],k)
# probability that k hps that are close enough form a consistent diagonal
def p_oo(k,p_hpSign):
	if k == 1:
		return 1
	else:
		return p_slash(k,p_hpSign) + p_backslash(k,p_hpSign)

# probability that in a window Wab of size la*lb there is at least one sb containing at least m hps (with no dispersed paralogy) spaced by gaps smaller or equal to g
# considered as the p-value of a sb
# Warning : an assumption is made that in the mhp there are only orthologies and no paralogies
def pValue(m,g,la,lb,nab,na,nb,p_hpSign):
	if min(nab,la,lb)+1 <= m: # means there are paralogies that make the computation impossible
		return None
	sum_ = 0
	for k in range(m,min(nab,la,lb)+1):
		sum__ = 0
		try:
			for i in range(m,k+1):
				sum__ += p_g_2D(i,g,la,lb)*p_oo(i,p_hpSign)
			sum_ += p_d(k,la,lb,nab,na,nb)*sum__
		except TypeError:
			sum_ = None
	return sum_

#from operator import mul    # or mul=lambda x,y:x*y
#from fractions import Fraction
#return float( int(reduce(mul, (Fraction(n-i, i+1) for i in range(p)), 1)) )
