# -*- coding: utf-8 -*-
# python 2.7 at least is needed
#
#                   PhylDiag v1.02
# 
# This code may be freely distributed and modified under the terms of the GNU General Public License version 3 or later and CeCILL. This should be distributed with the code. If you do not have a copy, see:
# 
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
# 
# Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team of the Institut de Biologie de l'Ecole Normale Supérieure and the individual authors. 
# 
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France
#
# This code is based on 
# ``PhylDiag : identifying complex cases of conserved synteny that include tandem duplications''
# authors :
# Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
# Matthieu MUFFATO (EBI, Cambridge, United Kingdom, muffato@ebi.ac.uk)
# Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)

# This code uses personal libraries in src/utils: myTools, myGenomes, myMaths, myFile, myPhylTree, mySvgDrawer.py, myProbas.py, myProteinTree and myDiags.py

##################
# PhylDiag v1.02 #
##################

###########
# INSTALL #
###########

From now on we consider that the user is in the root directory of the PhylDiag deposit 
: cd <PhylDiagFolder>

Install PhylDiag
-----------------
install python 2.7 :
: sudo apt-get install python2.7

Give excecution rights to scripts :
: chmod +x src/*.py
: chmod +x src/utils/*.py
-----------------

##############################################################
# preprocessing step : define gene families using gene trees #
##############################################################

Convert nhx (or .nwk, newick) gene trees to our tabular format (phylTree):
: src/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.phylTree

Convert a newick species tree into a phylTree species tree:
: src/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree

Extract the ancestral gene content (ancGene) from the gene trees:
: src/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.phylTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.phylTree

These ancGenes files can be used to define gene families
usually when two species Sa and Sb are compared, gene families are defined by ancGenes.<LCA(Sa,Sb)>

##########################
# processing of PhylDiag #
##########################

The core of PhylDiag is in src/utils/myDiags.py
Probability calculations are in src/utils/myProbas.py
The wrapper of PhylDiag is src/phylDiag.py

Access to the manual of PhylDiag :
: src/phylDiag.py

This returns :
//////////////////////
- ERROR - Not enough arguments
 Usage : src/phylDiag.py
	1: genome1 <type 'file'>
	2: genome2 <type 'file'>
	3: ancGenes <type 'file'>
	  -gapMax <type 'str'> (None)
	+/-sameStrand (True)
	  -filterType <type 'str'> (['None', 'InCommonAncestor', 'InBothSpecies'])
	  -minChromLength <type 'int'> (1)
	  -distanceMetric <type 'str'> (CD)
	  -pThreshold <type 'float'> (0.001)
	  -nbHpsRecommendedGap <type 'int'> (2)
	  -targetProbaRecommendedGap <type 'float'> (0.01)
	  -validateImpossToCalc_mThreshold <type 'int'> (3)
	+/-multiProcess (True)
	+/-verbose (False)


	Wrapper for PhylDiag Library
//////////////////////
numbered parameters are required,
other parameters are optional (the default value is between brackets)

3 files are required :
-2 genomes
-1 ancGene file (for gene family definitions)
Examples are given in <PhylDiagRootPath>/data: human and mouse genomes and the anGene defining families from the Euarchontoglire ancestor (LCA(human,mouse))

By default gapMax is set to 'None', thus PhylDiag chooses itself the advised gapMax (see article)
The distance metric may be either the 'DPD', 'ED', 'MD' or 'CD'; the default value is 'CD'
pThreshold is the p-value threshold for the statistical validation of synteny blocks
to follow step 1 of PhylDiag as described in the article, users should add -filterType=InCommonAncestor to filter lineage specific gene apparitions

Other parameters are not explained in PhylDiag article and users should avoid changing them 

A standard way to launch PhylDiag is thus:
: src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -filterType=InCommonAncestor > res/syntenyBlocks.txt
(this way, the gapMax is set with the advised value. Take care that the default distance metric is the chebyshev distance metric (CD))

A more customised way is for instance:
: src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -gapMax=11 -distanceMetric=MD -filterType=InCommonAncestor > res/syntenyBlocks_MD11.txt

Adding "+verbose" (set the verbose boolean to True) allows more informations in the error log, thus users willing to have more informations should launch PhylDiag like this:
: src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -distanceMetric=DPD -gapMax=5 -pThreshold=0.001 +verbose -filterType=InCommonAncestor > res/syntenyBlocks_DPD11.txt 2> res/logErr_DPD11.txt

Adding the "-multiProcess" (set the multiProcess boolean to False) option allow to launch phylDiag in sequential mode on one process alone. This may be useful for users that wants to launch phylDiag on a personal computer. Adding "nice -n 19" before the command line "phylDiag ..." may also be useful in order to prevent any crash of the OS if the CPU or memory consumption becomes greedy

"-nbHpsRecommendedGap=X" sets the number of homologies to X for the target synteny block in the average MHP for the recommended gapMax type (the default value is 2) and "-targetProbaRecommendedGap" is the value of the p-value of this target synteny block of "nbHpsRecommendedGap" hps. See the Supp Data of PhylDiag for more explicit informations

"-minChromLength" is the minimal length of considered chromosome for pairwise comparisons. Using a minChromLength of 2 may improve a lot the speed of the algorithm without loosing any sb (since sb comes from chromosomes of at least 2 genes) but may give wrong statistics on genomes at the beginning of the error log

###################
# PhylDiag Viewer #
###################

Access to the manual of PhylDiag's viewer :
: src/phylDiagHomologyMatrixViewer.py

This returns :
////////////////////////////////
- ERROR - Not enough arguments
 Usage : src/phylDiagHomologyMatrixViewer.py
	1: genome1 <type 'file'>
	2: genome2 <type 'file'>
	3: ancGenes <type 'file'>
	4: chr1:deb1-fin1 <type 'str'>
	5: chr2:deb2-fin2 <type 'str'>
	  -gapMax <type 'str'> (None)
	+/-consistentSwDType (True)
	  -filterType <type 'str'> (InCommonAncestor)
	  -minChromLength <type 'int'> (1)
	  -pThreshold <type 'float'> (0.001)
	  -out:SyntenyBlocks <type 'str'> (./res/syntenyBlocksDrawer.txt)
	+/-mode:chromosomesRewrittenInTbs (False)
	+/-convertGenicToTbCoordinates (False)
	  -distanceMetric <type 'str'> (CD)
	  -nbHpsRecommendedGap <type 'int'> (2)
	  -targetProbaRecommendedGap <type 'float'> (0.01)
	  -out:ImageName <type 'str'> (./res/homologyMatrix.svg)
	+/-verbose (True)


	Show the homology matrix with coloured synteny blocks (also called diagonals).
	- Each colour represents a synteny block.
	- On the x-axis are the genes of the 1st genome in the desired window
	- On the y-axis are the genes of the 2nd genome in the desired window
	- Each coloured rectangle in the matrix represents a filiation relationship restricted to the ancGene species. Couloured rectangle means the corresponding horizontal and vertical genes come from the same ancestral gene in the ancGene species. Take care that coulour genes are not homology relationship. For instance if the ancestor had two paralogs in its genome, there will be two distinct filiation relationships for genes herited from these genes. Genes are coloured if they are in the same filiation relationship class.
	- '+' indicates filiation relationships that have horizontal gene and vertical gene in the same direction on both chromosomes
	  '-' indicates filiation relationships that have horizontal gene and vertical gene in opposite directions on each chromosomes
///////////////////////////////////////////////////////////////////////////////////
numbered parameters are required
other parameters are optional  (the default value is between brackets)

To see the Matrix of Homologies (MH) of the human X chromosome compared to the mouse X chromosome, launch:
: src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -gapMax=5 -distanceMetric=MD -out:ImageName=res/MH_MD5.svg -out:SyntenyBlocks=res/syntenyBlocksDrawerMH_MD5.txt

"X:1-~" means X chromosome from the 1st gene to the last gene
For instance "4:45-80" means the 4th chromosome from the 45th gene to the 80th gene

The image MH.svg can be viewed with an internet browser as firefox 

It is also possible to draw the Matrix of Homology Packs (MHP):
: src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ +mode:chromosomesRewrittenInTbs -distanceMetric=MD -gapMax=5 -out:ImageName=./res/MHP_MD5.svg -out:SyntenyBlocks=./res/syntenyBlocksDrawerMHP_MD5.txt

Many parameters can be customised, for instance a user can run:
: Title=PhylDiag && S1=Homo.sapiens && S2=Mus.musculus && C1=X && R1="100-250" && C2=X && R2="1-100" && DM="DPD" && D=10 && src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${DM} -gapMax=${D} $C1:$R1 $C2:$R2 -out:ImageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_MHP.svg -out:SyntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001
