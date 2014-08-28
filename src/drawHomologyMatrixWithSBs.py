#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

"""
Build the homology matrix with synteny blocks using the mySvgDrawer Library
"""

import sys
import collections
import itertools
import random
import utils.myTools
import utils.mySvgDrawer as svgDrw

# draw either the mh or the mhp, if draw mode is 'writeinTB'
# inputs :
# 	genesStrandsCX = [+1, -1, ...] of length = to nX
# 	genesRemovedDuringFilteringC1 = [..., i, ...] with i the index of the gene removed during the filtering process (CX : Chromosome X)
# 	tbWithNoHomologyInWindowC1 = [..., [i6,i7,i8], ...] list of tbs with no homologies in the window, inside the index of genes. If draw mode is 'writeinTB' : tbWithNoHomologyInWindowC1 = [..., i6, ...] : just the index of the TB
# 	hpSigns = { i1 : {i2 : s, ...}...} with hpSign the hp sign (s1*s2) of the homology at the (i1,i2) coordinate (i1-th gene on C1 and i2-th gene on C2)
# 	homologyGroupsInWindow = [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...] a list of 2-uples of homologous tbs indices to be coloured in the same color
# 		with tb_CX_X : [..., i, ...] list of indices of genes in the tb
# 	diagIndices = [..., [...,(i16,j16),...], ...] list of diagonals with diagonals = list of all the points of the diagonal
# output :
#	string with the svg drawing of the mhp (or the mh)
def drawHomologyMatrix(((begC1,endC1),(begC2,endC2)), (genesStrandsC1, genesStrandsC2), (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2), (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2), hpSigns, homologyGroupsInWindow, diagsIndices, outputFileName=None, maxWidth=100, maxHeight=100 , symbolsInGenes=None):
	begC1=begC1+1 # For the print on screen
	begC2=begC2+1
	endC1=endC1+1
	endC2=endC2+1
	# nx = number of genes on the first genome
	# ny = number of genes on the second genome
	nx = len(genesStrandsC1)
	ny = len(genesStrandsC2)

	# corresponds to a palete of colors in the css file
	availableColors=range(8,40,2) # step = 2
	availableColors.reverse()
	availableGreys = range(0,15,1)
	random.shuffle(availableGreys)
	availableColorsForDiags = range(2,43)
	random.shuffle(availableColorsForDiags)

	# the size of the components of the matrix is chosen using the smallest and more restricting dimension (contains more genes comparing to its size)
	sizeCase = float (      min( float(maxWidth) / (nx + 3), float(maxHeight) / (ny + 3))    )# +3 for margins
	sizeText = float( sizeCase*0.9 )
	width = float(nx+3) * sizeCase
	height = float(ny+3) * sizeCase
	print >> sys.stderr, "width=", width
	print >> sys.stderr, "height=", height

	scene = svgDrw.Scene(name='homology_matrix', width=width, height=height)

	nbLinesX= nx+1 #Nb of vertical lines (x varies) in the matrix
	nbLinesY= ny+1 #Nb of horizontal lines (y varies) in the matrix

	# draw lines of chromosomes
	# offset_genes : corresponds to the chromosome lines positions passing through the middle of the genes
	offset_genes_x = sizeCase + sizeCase/2
	offset_genes_y = sizeCase + sizeCase/2
	scene.add(svgDrw.Line((offset_genes_x, height - offset_genes_y), ( width-sizeCase, height - offset_genes_y), width=0.1*sizeCase))
	scene.add(svgDrw.Line((offset_genes_x, sizeCase), ( offset_genes_x, height - (offset_genes_y)), width=0.1*sizeCase))
	offset_matrix_x = 2*sizeCase
	offset_matrix_y = 2*sizeCase
	# draw Diagonals first because they are on the background
	print >> sys.stderr, "Nb of diagonals showed = ", len(diagsIndices)
	for diag in diagsIndices:
		#color = availableGreys.pop() if len(availableGreys) > 0 else random.randint(0,230)
		color = availableColorsForDiags.pop() if len(availableColorsForDiags)>0 else random.choice(range(2,43))
		for (i,j) in diag:
			cx_s=i*sizeCase
			cy_s=j*sizeCase
			scene.add(svgDrw.Rectangle((cx_s+offset_matrix_x,height-(cy_s+sizeCase+offset_matrix_y)), sizeCase, sizeCase, fill_opacity=0.90, svgClass = "HomologGroup%s" % color ))
			#scene.add(scgDrw.Rectangle((cx_s,height-(cy_s+sizeCase)), sizeCase, sizeCase, fill_opacity=0.90, svgClass = "chromgrp%s" % color ))
			#scene.add(svgDrw.Rectangle((cx_s,height-(cy_s+sizeCase)), sizeCase, sizeCase, fill_opacity=fill_opacity=0.90, color=(color,color,color))
		# draw rectangles around diagonals
		min_i = min(diag,key=lambda x:x[0])[0]
		max_i = max(diag,key=lambda x:x[0])[0]
		min_j = min(diag,key=lambda x:x[1])[1]
		max_j = max(diag,key=lambda x:x[1])[1]
		cx_s=min_i*sizeCase
		cy_s=max_j*sizeCase
		scene.add(svgDrw.Rectangle((cx_s+offset_matrix_x,height-(cy_s+offset_matrix_y+sizeCase)), (max_j-min_j)*sizeCase + sizeCase, (max_i-min_i)*sizeCase + sizeCase,  stroke='black', fill='none', strokeWidth=0.2*sizeCase ))

	# tick lines
	widthTicks = min(float(width)/1000, float(height)/1000)
	sizeTextTicks = widthTicks*10
	for i,ni in enumerate(range(begC1,endC1)): # TODO : better place ticks
		cx = i*sizeCase
		if ni%10==0:
			scene.add(svgDrw.Line((offset_matrix_x+sizeCase/2+cx, height - offset_genes_y/2), (offset_matrix_x+sizeCase/2+cx, height - offset_genes_y), width=widthTicks))
		if ni%50 == 0:
			cxText = offset_matrix_x+sizeCase/2+cx
			cyText = height - max(offset_genes_y/2, sizeTextTicks/2)
			cxx = offset_matrix_x+cx+sizeCase/2
			if nx > 750 or ny > 750:
				if ni%100 == 0:
					scene.add(svgDrw.Text((cxText, cyText), str(ni), text_anchor="middle", size=sizeTextTicks))
					scene.add(svgDrw.Line((cxx, height - offset_matrix_y), (cxx, sizeCase), width=sizeCase*0.1))
			else:
				scene.add(svgDrw.Text((cxText, cyText), str(ni), text_anchor="middle", size=sizeTextTicks))
				scene.add(svgDrw.Line((cxx, height - offset_matrix_y), (cxx, sizeCase), width=sizeCase*0.1))

	for j,nj in enumerate(range(begC2,endC2)): # TODO : better place ticks
		cy = j*sizeCase
		if nj%10==0:
			scene.add(svgDrw.Line((offset_genes_x/2, height - (offset_matrix_y+sizeCase/2+cy)), (offset_genes_x, (height - (offset_matrix_y+sizeCase/2+cy))), width=widthTicks))
		if nj%50 == 0:
			cyText = height - (offset_matrix_y+sizeCase/2+cy)
			cxText = max(offset_genes_x/2, sizeTextTicks/2)
			cxx = offset_genes_x/2
			cyy = height - (offset_matrix_y+sizeCase/2+cy)
			if nx > 750 or ny > 750:
				if nj%100 == 0:
					#scene.add(svgDrw.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCase,cxx+sizeCase,cyy)))
					scene.add(svgDrw.Text((cxText, cyText), str(nj), text_anchor="middle", size=sizeTextTicks, transform="rotate(90,%s,%s)" % (cxText,cyText)))
					scene.add(svgDrw.Line((offset_matrix_x, cyy), (width-sizeCase, cyy), width=sizeCase*0.1))
			else:
				#scene.add(svgDrw.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCase,cxx+sizeCase,cyy)))
				scene.add(svgDrw.Text((cxText, cyText), str(nj), text_anchor="middle", size=sizeTextTicks, transform="rotate(90,%s,%s)" % (cxText,cyText)))
				scene.add(svgDrw.Line((offset_matrix_x, cyy), (width-sizeCase, cyy), width=sizeCase*0.1))

	if nx < 300 and ny < 300:
		for i in range(nbLinesX):
			cxLine = i*sizeCase
			scene.add(svgDrw.Line((cxLine+offset_matrix_x, height - offset_matrix_y), (cxLine+offset_matrix_x, height - ((nbLinesY-1)*sizeCase + offset_matrix_y)), width=sizeCase*0.01))
		for j in range(nbLinesY):
			cyLine = j*sizeCase
			scene.add(svgDrw.Line((offset_matrix_x, height - (cyLine+offset_matrix_y) ), (offset_matrix_x + (nbLinesX-1)*sizeCase, height - (cyLine+offset_matrix_y) ), width=sizeCase*0.01))

		chromosome1={}
		chromosome2={}

		# create chromosomes
		for i1,s1 in enumerate(genesStrandsC1):
			cx = i1*sizeCase
			symbol = str(symbolsInGenes[0][i1]) if symbolsInGenes != None else None
			chromosome1[i1]=svgDrw.Gene((cx+offset_matrix_x, height-offset_genes_y), (cx+sizeCase+offset_matrix_x, height-offset_genes_y), strand=s1, width=sizeCase*0.7, stroke_width=0.05*sizeCase, SVGclass=None, text=symbol)
		for i2,s2 in enumerate(genesStrandsC2):
			cy = i2*sizeCase
			symbol = str(symbolsInGenes[1][i2]) if symbolsInGenes != None else None
			chromosome2[i2]=svgDrw.Gene((offset_genes_x,height-(cy+offset_matrix_y)), (offset_genes_x, height-(cy+sizeCase+offset_matrix_y)), strand=s2, width=sizeCase*0.7, stroke_width=0.05*sizeCase, SVGclass=None, text=symbol)

		# give a color to each gene using homology relationships
		for (tbs1,tbs2) in homologyGroupsInWindow:
			color = availableColors.pop() if len(availableColors) > 0 else random.choice(range(8,40,2))
			for tb1 in tbs1:
				if isinstance(tb1,list):
					for i1 in tb1:
						chromosome1[i1].SVGclass = "HomologGroup%s" % color
				else:
					chromosome1[tb1].SVGclass = "HomologGroup%s" % color
			for tb2 in tbs2:
				if isinstance(tb2,list):
					for i2 in tb2:
						chromosome2[i2].SVGclass = "HomologGroup%s" % color
				else:
					chromosome2[tb2].SVGclass = "HomologGroup%s" % color

		# give grey levels to genes that have no homology in the window
		for tb in tbWithNoHomologyInWindowC1:
			grey = availableGreys.pop() if len(availableGreys) > 0 else random.choice(range(0,13,1))
			if isinstance(tb,list):
				for i1 in tb:
					chromosome1[i1].SVGclass = "NoHomologyInWindow%s" % grey
			else:
				chromosome1[tb].SVGclass = "NoHomologyInWindow%s" % grey
		for tb in tbWithNoHomologyInWindowC2:
			grey = availableGreys.pop() if len(availableGreys) > 0 else random.choice(range(0,13,1))
			if isinstance(tb,list):
				for i2 in tb:
					chromosome2[i2].SVGclass = "NoHomologyInWindow%s" % grey
			else:
				chromosome2[tb].SVGclass = "NoHomologyInWindow%s" % grey

		for i1 in genesRemovedDuringFilteringC1:
			chromosome1[i1].SVGclass = "SpeciesSpecificGenes"
		for i2 in genesRemovedDuringFilteringC2:
			chromosome2[i2].SVGclass = "SpeciesSpecificGenes"

		for i1 in chromosome1:
			scene.add(chromosome1[i1])
		for i2 in chromosome2:
			scene.add(chromosome2[i2])

		# fill homologies with +1,-1 or ? or 0
		nonZeroValues=[]
		for i1 in hpSigns:
			for i2 in hpSigns[i1]:
				nonZeroValues.append((i1,i2))
				#s = hpSigns[i1][i2][0][1]
				s = hpSigns[i1][i2]
				cx = i1*sizeCase + float(sizeCase)/2
				cy = i2*sizeCase + float(sizeCase)/2
				assert s == +1 or s== -1 or s == None, "s=%s" % s
				assocValue = (("+" if s == +1 else "-") if s != None else '?')
				scene.add(svgDrw.Text((cx+offset_matrix_x,height-(cy+sizeText*0.16+offset_matrix_y)), assocValue, text_anchor="middle", size=sizeText))
		if nx < 20 and ny < 20:
			for (i1,i2) in itertools.product(range(nx),range(ny)):
				if (i1,i2) not in nonZeroValues:
					cx = i1*sizeCase + float(sizeCase)/2
					cy = i2*sizeCase + float(sizeCase)/2
					scene.add(svgDrw.Text((cx+offset_matrix_x,height-(cy+sizeText*0.16+offset_matrix_y)), "0", text_anchor="middle", size=sizeText, fill=(200,200,200), stroke=None))

	else:
		# represent homologies with a black rectangle
		for i1 in hpSigns:
			for i2 in hpSigns[i1]:
				if (i1,i2) not in [dot for diag in diagsIndices for dot in diag]:
					cx_s=i1*sizeCase
					cy_s=i2*sizeCase
					scene.add(svgDrw.Rectangle((cx_s+offset_matrix_x,height-(cy_s+sizeCase+offset_matrix_y)), sizeCase, sizeCase, fill=(0,0,0), fill_opacity=0.90))
		print >> sys.stderr, "Warning : some supplementary informations are not displayed because one of the two dimension of the window is > 300"


	if outputFileName != None:
	        scene.write_svg(filename=str(outputFileName))
	#scene.display()
	return scene.strarray()

def test(outFileName):
    scenario = arguments["scenario"]
    if scenario==1 :
	nx=11
    	ny=7
    	homologySameStrand=[(1,1), (2,2), (3,3), (4,4), (5,5), (8,3), (9,4)]
    	homologyOppositeStrand=[]
	diags=[[(1,1), (2,2), (3,3), (4,4), (5,5)],[(8,3), (9,4)]]
    elif scenario ==2:
    	nx=9
    	ny=9
    	homologySameStrand=[(1,1), (2,2), (3,4), (3,5), (4,3), (5,3), (6,6), (7,7)] # not necessary to mark all for hps
    	homologyOppositeStrand=[(3,3), (4,4), (4,5), (5,4), (5,5)]
    	diags=[homologySameStrand + homologyOppositeStrand]
    elif scenario == 3:
      	nx=6
    	ny=6
    	homologySameStrand=[(1,1), (2,2)]
    	homologyOppositeStrand=[(3,3), (4,4)]
    	diags=[homologySameStrand, [(3,3)], [(4,4)]]
    elif scenario == 4:
    	nx=10
    	ny=7
    	homologySameStrand=[(3,3)]
    	homologyOppositeStrand=[(1,5), (2,4), (4,3), (5,3), (6,3), (7, 2), (8, 1)]
    	diags=[homologySameStrand + homologyOppositeStrand]
    elif scenario == 5:
    	nx=9
    	ny=9
    	homologySameStrand=[(1,1), (2,2), (3,3), (6,6), (7,7)]
    	homologyOppositeStrand=[]
    	diags=[homologySameStrand, homologyOppositeStrand]
    elif scenario == 6:
    	nx=13
    	ny=13
    	homologySameStrand=[(1,6), (2,7), (4,8), (5,9), (6,10), (10,8), (11,9)]
    	homologyOppositeStrand=[(7,3),(8,2), (3,8)]
    elif scenario == 7:
    	# For the paper on the method to extract diagonals
    	nx=8
    	ny=9
    	homologySameStrand=[(1,2), (2,1), (3,3), (4,4), (4,5), (4,6), (5,4), (5,5), (5,6), (6,7)]
    	homologyOppositeStrand=[(1,1),(2,2)]
    	diags=[]
    elif scenario == 8:
    	# For the paper on the method to extract diagonals
    	nx=6
    	ny=6
    	homologySameStrand=[(1,1), (2,2), (3,3), (4,4)]
    	homologyOppositeStrand=[]
    	diags=[]
    elif scenario == 9:
	# For the paper on the method to extract diagonals
    	nx=7
    	ny=5
    	homologySameStrand=[(2,2)]
    	homologyOppositeStrand=[(1,3),(3,2),(4,2),(5,1)]
    	diags=[]
    elif scenario == 10:
    	# For the paper on the method to extract diagonals
    	nx=5
    	ny=5
    	homologySameStrand=[]
    	homologyOppositeStrand=[(1,3),(2,2),(3,1)]
    	diags=[]
    elif scenario == 11:
    	# For the paper on the method to extract diagonals (merge and probabilities)
    	nx=11
    	ny=8
	homologyPlusAssociatedValues=[(1,(5,))]
    	homologyMinusAssociatedValues=[(4,(4,)),(5,(3,)),(9,(1,))]
	homologyUnknownAssociatedValues=[(7,(6,)),(8,(2,))]
	hpSigns = collections.defaultdict(dict)
	genesStrandsC1=[]
	genesStrandsC2=[]
	for i1 in range(nx):
		genesStrandsC1.append(random.choice([-1,1]))
	for i2 in range(ny):
		genesStrandsC2.append(random.choice([-1,1]))

	for (indexesAVs,aV) in itertools.izip([homologyPlusAssociatedValues, homologyMinusAssociatedValues, homologyUnknownAssociatedValues],[+1,-1,None]):
		for (i1,i2s) in indexesAVs:
			for i2 in i2s:
				hpSigns[i1][i2]=aV
				# Determine possible orientations
				strand = random.choice([-1,+1])
				if  aV == +1:
					genesStrandsC1[i1]=strand
					genesStrandsC2[i2]=strand
				elif aV == -1:
					genesStrandsC1[i1]=strand
					genesStrandsC2[i2]=-strand
				elif aV == None:
					b=random.choice([True, False])
					genesStrandsC1[i1]=None if b else strand
					genesStrandsC2[i2]=strand if b else None
				else:
					raise ValueError()
	begC1,endC1 = 0,nx-1
	begC2,endC2 = 0,ny-1
	tbWithNoHomologyInWindowC1=[0,2,3,6,10]
	tbWithNoHomologyInWindowC2=[0,7]
	homologyGroupsInWindow=[([1],[5]),([4],[4]),([5],[3]),([9],[1]),([7],[6]),([8],[2])]
    	diagsIndices=[[(4,4),(5,3),(8,2),(9,1)]]
	genesRemovedDuringFilteringC1=[]
	genesRemovedDuringFilteringC2=[]
	symbolsInGenes=[[1,1,1,1,1,1,1,2,2,1,1],[1,1,2,1,1,1,2,1]]

	drawHomologyMatrix(((begC1,endC1),(begC2,endC2)), (genesStrandsC1, genesStrandsC2), (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2), (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2), hpSigns, homologyGroupsInWindow, diagsIndices, outputFileName=outFileName, maxWidth=100, maxHeight=100 , symbolsInGenes=symbolsInGenes)

if __name__ == '__main__':
	arguments = utils.myTools.checkArgs([("scenario",int)],[("out:FileName",str,"image.svg")],__doc__)
	test(arguments['out:FileName'])

