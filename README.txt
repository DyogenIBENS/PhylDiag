 ____  _           _ ____  _
|  _ \| |__  _   _| |  _ \(_) __ _  __ _
| |_) | '_ \| | | | | | | | |/ _` |/ _` |
|  __/| | | | |_| | | |_| | | (_| | (_| |
|_|   |_| |_|\__, |_|____/|_|\__,_|\__, |
             |___/                 |___/

PhylDiag version 2.01 (02/06/2017)
python v2.7 at least is needed

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:
    1) LICENSE-GPL.txt (or http://www.gnu.org/licenses/gpl-3.0-standalone.html)
    2) LICENCE-CeCILL.txt (or http://www.cecill.info/licences/Licence_CeCILL_V2-en.html)

Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure (IBENS) 46 rue d'Ulm Paris, and the individual authors.

Copyright © 2015 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
mail : hrc@ens.fr or jlucas@ens.fr

This code is based on two publications
[1] "PhylDiag : identifying complex synteny blocks that include tandem duplications using phylogenetic gene trees"
    Authors :
    Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
    Matthieu MUFFATO (EBI, Cambridge, United Kingdom, muffato@ebi.ac.uk)
    Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)
    Published in BMC Bioinformatics, August 2014
FIXME update as soon as published
[2] "High precision detection of conserved segments from synteny blocks"
    Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
    Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)
    submitted in PLOS ONE, minor revisions

Both articles are in the doc/ folder, and they can be downloaded at the next urls:
Pdf version:
http://www.biomedcentral.com/content/pdf/1471-2105-15-268.pdf
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155083/pdf/12859_2014_Article_6554.pdf
FIXME 2nd article links
Online article:
http://www.ncbi.nlm.nih.gov/pubmed/?term=phyldiag
http://www.biomedcentral.com/1471-2105/15/268
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155083/

Supplementary data are also in the doc/ folder or can be accessed at the next url:
http://www.biomedcentral.com/content/supplementary/1471-2105-15-268-S1.pdf

This code uses a personal library: contained in git@depot.biologie.ens.fr:LibsDyogen

###########
###########
# INSTALL #
###########
###########

# Install the LibsDyogen library first
# Read : https://raw.githubusercontent.com/DyogenIBENS/LibsDyogen/master/README.txt
# From now on we assume that the path to the folder libsDyogen is in the PYTHONPATH

# The easiest way to install PhylDiag is to launch the remote script INSTALL.sh hosted on github
# This script will clone the github deposit itself
# The installation will be set in /home/${USER}/Libs/PhylDiag

# Install curl, if you don't have it
sudo apt-get update
sudo apt-get install curl
# Use curl to execute the remote file INSTALL.sh hosted on github
bash <(curl -s https://raw.githubusercontent.com/DyogenIBENS/PhylDiag/master/INSTALL.sh)

# If this does not work, follow the next indications

# Install dependencies:
# LibsDyogen

# Choose a path for the parent folder of PhylDiag (here it is /home/<user>/Libs)
PATH_PARENT_PHYLDIAG="/home/${USER}/Libs"
mkdir -p ${PATH_PARENT_PHYLDIAG}
cd ${PATH_PARENT_PHYLDIAG}
PATH_PHYLDIAG=${PATH_PARENT_PHYLDIAG}/PhylDiag
# Clone the PhylDiag deposit
git clone https://github.com/DyogenIBENS/PhylDiag ${PATH_PHYLDIAG}
# If necessary give execution rights
chmod +x ${PATH_PHYLDIAG}/src/*.py
chmod +x ${PATH_PHYLDIAG}/src/analysis/*.py
chmod +x ${PATH_PHYLDIAG}/src/postprocessing/*.py

# It should be installed. Check that everything works properly by launching :
cd PhylDiag
bash ./checkPhylDiagIntegrity.sh
# this should run several tests on phyldiag

##########
##########
# MANUAL #
##########
##########

##############################################################
# Preprocessing steps: define gene families using gene trees #
##############################################################

# the scripts for pre-processing species tree, gene trees, gene families... have been moved to LibsDyogen/scripts
# because they are used by other programs than PhylDiag

##########################
# Processing of PhylDiag #
##########################

The core of PhylDiag is in LibsDyogen/utils/myDiags.py
Probability calculations are in LibsDyogen/utils/myProbas.py
The wrapper of PhylDiag is src/phylDiag.py

# Access to the synthesis of the manual of PhylDiag command line:
src/phylDiag.py

This returns :
///////////////////////////////////////////////////////////////////////////////
- ERROR - Not enough arguments
 Usage : src/phylDiag.py
	1: genome1 <type 'file'>
	2: genome2 <type 'file'>
	3: families <type 'file'>
	  -filterType <type 'str'> (InBothGenomes)
	  -tandemGapMax <type 'int'> (10)
	  -gapMax <type 'str'> (5)
	+/-distinguishMonoGenicDiags (True)
	  -distanceMetric <type 'str'> (CD)
	  -pThreshold <type 'str'> (None)
	  -gapMaxMicroInv <type 'str'> (1)
	+/-identifyMonoGenicInvs (True)
	+/-identifyMicroRearrangements (True)
	  -truncationMax <type 'str'> (10)
	  -minChromLength <type 'int'> (2)
	+/-sameStrand (True)
	  -nbHpsRecommendedGap <type 'int'> (2)
	  -targetProbaRecommendedGap <type 'float'> (0.01)
	  -validateImpossToCalc_mThreshold <type 'int'> (3)
	  -optimisation <type 'str'> (cython)
	+/-verbose (False)
	+/-removeUnofficialChromosomes (True)


Wrapper for PhylDiag Library
///////////////////////////////////////////////////////////////////////////////
Numbered parameters are required and other parameters are optional, used default value are written between brackets.

3 files are required :
-2 genomes
-1 ancGenes file (for gene family definitions)
Examples of genomes and ancGenes files are given in <PhylDiagRootPath>/data:
human and mouse genomes and the anGene defining families from the Euarchontoglire ancestor = MRCA(human,mouse)

FIXME explain all options
* By default gapMax is set to 'None', thus PhylDiag chooses itself the advised gapMax, see article [1].
* The distance metric may be either the 'DPD', 'ED', 'MD' or 'CD'; the default value is 'CD'
* pThreshold is the p-value threshold for the statistical validation of synteny blocks, default value is None, meanning
that their is no statistical validation of synteny blocks.
* By default the filtering of extant genomes is InBothGenomes, meanning that only homologs are kept.

Other parameters are not explained in PhylDiag article and novice users should avoid changing them

# A standard way to launch PhylDiag is thus:
src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -tandemGapMax=5 -gapMax=5 -truncationMax=5 > res/syntenyBlocks.txt

# Adding "+verbose" (set the verbose boolean to True) returns more information in the logErr:
src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -tandemGapMax=5 -gapMax=5 -truncationMax=5 +verbose > res/syntenyBlocks.txt 2> res/logErr.txt

###################
# PhylDiag Viewer #
###################

# Access to the synthesis of the manual of PhylDiag command line:
src/phylDiagHomologyMatrixViewer.py

This returns :
/////////////////////////////////////////////////////////////////////////
- ERROR - Not enough arguments
 Usage : src/phylDiagHomologyMatrixViewer.py
	1: genome1 <type 'file'>
	2: genome2 <type 'file'>
	3: ancGenes <type 'file'>
	4: chr1:deb1-fin1 <type 'str'>
	5: chr2:deb2-fin2 <type 'str'>
	  -gapMax <type 'str'> (None)
	+/-sameStrand (True)
	  -filterType <type 'str'> (InCommonAncestor)
	  -minChromLength <type 'int'> (1)
	  -pThreshold <type 'float'> (0.001)
	  -out:syntenyBlocks <type 'str'> (./res/syntenyBlocksDrawer.txt)
	+/-mode:chromosomesRewrittenInTbs (False)
	+/-convertGenicToTbCoordinates (False)
	  -distanceMetric <type 'str'> (CD)
	  -nbHpsRecommendedGap <type 'int'> (2)
	  -targetProbaRecommendedGap <type 'float'> (0.01)
	  -out:imageName <type 'str'> (./res/homologyMatrix.svg)
	+/-verbose (True)


	Show the homology matrix with coloured synteny blocks (also called diagonals).
	- Each colour represents a synteny block.
	- On the x-axis are the genes of the 1st genome in the desired window
	- On the y-axis are the genes of the 2nd genome in the desired window
	- Each coloured rectangle in the matrix represents a filiation relationship restricted to the ancGene species. Couloured rectangle means the corresponding horizontal and vertical genes come from the same ancestral gene in the ancGene species. Take care that coulour genes are not homology relationship. For instance if the ancestor had two paralogs in its genome, there will be two distinct filiation relationships for genes herited from these genes. Genes are coloured if they are in the same filiation relationship class.
	- '+' indicates filiation relationships that have horizontal gene and vertical gene in the same direction on both chromosomes
	  '-' indicates filiation relationships that have horizontal gene and vertical gene in opposite directions on each chromosomes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
Numbered parameters are required and other parameters are optional, used default value are written between brackets.

# To see the Matrix of Homologies (MH) of the human X chromosome compared to the mouse X chromosome, launch:
src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -tandemGapMax=5 -gapMax=5 -truncationMax=5 -out:imageName=res/MH.svg -out:syntenyBlocks=res/syntenyBlocksDrawerMH.txt

"X:1-~" means X chromosome from the 1st gene to the last gene
For instance "4:45-80" means the 4th chromosome from the 45th gene to the 80th gene

The image MH.svg can be viewed with an internet browser as firefox.

# It is also possible to draw the Matrix of Homology Packs (MHP):
src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ +mode:chromosomesRewrittenInTbs -tandemGapMax=5 -gapMax=5 -truncationMax=5 -out:imageName=./res/MHP.svg -out:syntenyBlocks=./res/syntenyBlocksDrawerMHP.txt

# Many parameters can be customised, for instance a user can run:
Title=PhylDiag && S1=Homo.sapiens && S2=Mus.musculus && C1=X && R1="100-250" && C2=X && R2="1-100" && DM="DPD" && D=10 && src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${DM} -gapMax=${D} $C1:$R1 $C2:$R2 -out:imageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_MHP.svg -out:syntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001

