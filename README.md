# PhylDiag
[![DOI](https://zenodo.org/badge/19742670.svg)](https://zenodo.org/badge/latestdoi/19742670)


From the comparison of two extant genomes and corresponding gene trees (or gene families), PhylDiag detects conserved segments, i.e. segments of chromosomes unbroken during evolution.

**Inputs**
1. two extant *genomes*, G1 and G2
2. *gene families*, F

A genome is a set of chromosomes.
A chromosome is a list of genes.
A gene is a pair (gene name, strand).

F is an associative array that links, for each family,
* the family name (**key**)
* to the set of names of the descendant genes (**values**)

The name of the family is often the name of the ancestral gene, at the root of the gene family.

If you use phylogenetic gene trees, utils in LibsDyogen/scripts  can
convert your trees into gene families.

**Outputs**
* Either (depending on the options of PhylDiag)
    * *conserved segments*
    * or *synteny blocks*, if
        * the identification of micro-rearrangements,
        * identification of mono-genic conserved segments,
        * identification of mono-genic inversions

        are disabled

PhylDiag is explained in more details in two publications
1. [PhylDiag : identifying complex synteny blocks that include tandem duplications using phylogenetic gene trees](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-268)
2. High precision detection of conserved segments from synteny blocks (FIXME update as soon as published)

and in a french thesis manuscript

3. [Étude de l’évolution de l’ordre des gènes de vertébrés par simulation](https://tel.archives-ouvertes.fr/tel-01398369/document)

##### Link between synteny blocks and conserved segments

> "there is no difference between a synteny block with no gap (g=0) and a conserved segment"
> -- <cite> [in High precision detection of conserved segments from synteny blocks]<cite/>

Remark: "no gap", corresponds to "no micro-rearrangement" in its context

*Conserved segments* can be considered as a specific type of *synteny blocks*.
For this reason you may see some *conserved segments* being named more generally *synteny blocks* in the code.

## Installation
[Install the LibsDyogen library first.](https://github.com/DyogenIBENS/LibsDyogen)
From now on we assume that the path to the folder LibsDyogen is in the PYTHONPATH.

The easiest way to install PhylDiag is to launch the remote script INSTALL.sh hosted on github.
This script will clone the github deposit itself.
The installation will be set in /home/${USER}/Libs/PhylDiag.

Install curl, if you don't have it
```
sudo apt-get update
sudo apt-get install curl
```
Use curl to execute the remote file INSTALL.sh hosted on github
```
bash <(curl -s https://raw.githubusercontent.com/DyogenIBENS/PhylDiag/master/INSTALL.sh)
```

If it did not work, follow the next instructions.

### Dependencies
Core dependencies
* LibsDyogen

### Detailed installation guidelines
Choose a path for the parent folder of PhylDiag (here it is /home/<user>/Libs)
```PATH_PARENT_PHYLDIAG="/home/${USER}/Libs"
mkdir -p ${PATH_PARENT_PHYLDIAG}
cd ${PATH_PARENT_PHYLDIAG}
PATH_PHYLDIAG=${PATH_PARENT_PHYLDIAG}/PhylDiag
```
Clone the PhylDiag deposit
```
git clone https://github.com/DyogenIBENS/PhylDiag ${PATH_PHYLDIAG}
```
If necessary give execution rights
```
chmod +x ${PATH_PHYLDIAG}/src/*.py
chmod +x ${PATH_PHYLDIAG}/src/analysis/*.py
chmod +x ${PATH_PHYLDIAG}/src/postprocessing/*.py
```

It should be installed.
You can verify that everything works properly with some tests
```
cd PhylDiag
bash ./checkPhylDiagIntegrity.sh
```

## Usage
### PhylDiag
* The core of PhylDiag is in LibsDyogen/utils/myDiags.py.
* Probability calculations are in LibsDyogen/utils/myProbas.py.
* The wrapper of PhylDiag is src/phylDiag.py.

Executing phylDiag without argument will show how to use the executable
```
src/phylDiag.py
```

It shows
```
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
```
Numbered parameters are required and other parameters are optional, default values used are written between brackets.

3 files are required :
- 2 genomes files
- 1 ancGenes file (for gene family definitions)

Examples of genomes and ancGenes files are given in <PhylDiagRootPath>/data.
For instance you could use human and mouse genomes and the anGenes defining families from the Euarchontoglire ancestor = MRCA(human,mouse)

(TODO explain all options)

* If the gapMax is set to 'None', PhylDiag chooses itself the advised gapMax, see article [1].
* The distance metric may be either the 'DPD', 'ED', 'MD' or 'CD'; the default value is 'CD'
* pThreshold is the p-value threshold for the statistical validation of synteny blocks, default value is 'None', meaning
that their is no statistical validation of synteny blocks.
* By default the filtering of extant genomes is InBothGenomes, meaning that only homologs are kept.

Other parameters are not explained in PhylDiag article and novice users should avoid changing them

A standard way to launch PhylDiag is thus
```
src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -tandemGapMax=5 -gapMax=5 -truncationMax=5 > res/consevedSegments.txt
```

Adding '+verbose' (set the verbose boolean to True) returns more information in the logErr
```
src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -tandemGapMax=5 -gapMax=5 -truncationMax=5 +verbose > res/conservedSegments.txt 2> res/logErr.txt
```

### PhylDiag Viewer

Here again, executing phylDiagHomologyMatrixViewer without argument will show how to use the executable
```
src/phylDiagHomologyMatrixViewer.py
```

This returns
```
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
```
Numbered parameters are required and other parameters are optional, used default value are written between brackets.

To see the Matrix of Homologies (MH) of the human X chromosome compared to the mouse X chromosome, execute
```
src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -tandemGapMax=5 -gapMax=5 -truncationMax=5 -out:imageName=res/MH.svg -out:syntenyBlocks=res/syntenyBlocksDrawerMH.txt
```
"X:1-~" means X chromosome from the 1st gene to the last gene
```
For instance "4:45-80" means the 4th chromosome from the 45th gene to the 80th gene
```
The image MH.svg can be viewed with an internet browser as firefox.

It is also possible to draw the Matrix of Homology Packs (MHP)
```
src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ +mode:chromosomesRewrittenInTbs -tandemGapMax=5 -gapMax=5 -truncationMax=5 -out:imageName=./res/MHP.svg -out:syntenyBlocks=./res/syntenyBlocksDrawerMHP.txt
```

Many parameters can be customised, for instance a user can execute
```
Title=PhylDiag && S1=Homo.sapiens && S2=Mus.musculus && C1=X && R1="100-250" && C2=X && R2="1-100" && DM="DPD" && D=10 && src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${DM} -gapMax=${D} $C1:$R1 $C2:$R2 -out:imageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_MHP.svg -out:syntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001
```

## Update
If you want to keep PhylDiag up to date, you need to update LibsDyogen first (see [the Update section of LibsDyogen](https://raw.githubusercontent.com/DyogenIBENS/LibsDyogen/master/README.md)).

Then
```
cd ${PATH_LIBSDYOGEN}
git pull
```
This will upgrade your local git deposit to the last commit.

If you want a more stable version, after `git pull`, you can downgrade to the latest tagged version (=stable release)
1. Get tags from the github deposit `git fetch --tags`
2. Get the latest tag name ``latestTag=$(git describe --tags `git rev-list --tags --max-count=1`)``
3. Checkout the latest tag `git checkout $latestTag`

Otherwise, after `git fetch --tags`
1. List all tagged versions: `git tag -l`
2. Checkout to the version you want: `git checkout <tagName>`


**Please ensure that the versions of PylDiag and LibsDyogen share the same tagged version or correspond to the last commits.**

## Contributing
If you want to contribute to this deposit please
1. fork it
2. create your feature branch: `git checkout -b my-new-feature`
3. commit your changes: `git commit -am 'Add some feature'`
4. push to the branch: `git push origin my-new-feature`
5. submit a pull request

## Credits
* Joseph Lucas: conceptualization and implementation of phylDiag
* Hugues Roest Crollius: supervision
* Lucas Tittmann: improved the clustering of tandem duplicates
* Nga thi thuy Nguyen: optimisation of the core algorithm of PhylDiag with cython
* Matthieu Muffato: implemented several python functions

## License
This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:
* LICENSE-GPL.txt (or http://www.gnu.org/licenses/gpl-3.0-standalone.html)
* LICENCE-CeCILL.txt (or http://www.cecill.info/licences/Licence_CeCILL_V2-en.html)
Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure (IBENS) 46 rue d'Ulm Paris, and the individual authors.

## Contacts

* [Joseph Lucas](jlucas@ens.fr)
* [Hugues Roest Crollius](hrc@ens.fr)
