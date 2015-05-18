#!/bin/bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

#############################################
#	Check integrity of pre-processing scripts #
#############################################

#preProcessCommandLines=(
## convet a .nhx tree into a protTree (forest of gene trees)
#"src/preprocessing/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.protTree"
## convet a .nwk tree into a phylTree
#"src/preprocessing/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree"
## extract ancGenes (family)  from the species tree and the forest of gene trees
#"src/preprocessing/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.protTree"
#)
#for line in "${preProcessCommandLines[@]}"
#	do
#		echo -e "${green}${line}${NC}"
#		eval ${line}
#done

#############################################
#	Check integrity of phylDiag.py					  #
#############################################

sbFileName(){
	if (( $# != 8 )); then
		    echo "Illegal number of parameters"
	fi
	s1=$1
	s2=$2
	a=$3
	# tandem gap max
	tgm=$4
	# gap max
	gm=$5
	# distance metric
	dm=$6
	#echo "${s1}_${s2}_f${a}_Tgm${tgm}Gm${gm}${dm}.sbs"
	# identify breakpoints within gaps
	ibwg=$7
	# overlap max
	om=$8
	if [ "${ibwg}" == "+" ]; then
		echo "${s1}_${s2}_f${a}_Tgm${tgm}Gm${gm}${dm}IbwgOm${om}.sbs"
	else
		echo "${s1}_${s2}_f${a}_Tgm${tgm}Gm${gm}${dm}Om${om}.sbs"
	fi
}

S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
phylDiagCommandLines=(
# phylDiag with default options
"src/phylDiag.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 > res/${S1}_${S2}_f${A}.sbs"
)
gm=5
dm=CD
phylDiagCommandLines+=(
# phylDiag with some options
"src/phylDiag.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 -filterType=InFamilies -distanceMetric=${dm} -gapMax=${gm} > res/${S1}_${S2}_f${A}_Gm${gm}${dm}.sbs"
"src/phylDiag.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 -filterType=InBothGenomes -distanceMetric=${dm} -gapMax=${gm} > res/${S1}_${S2}_fIBS${A}_Gm${gm}${dm}.sbs"
)
tgm=9
gm=5
dm=CD
ibwg='+'
om=10
phylDiagCommandLines+=(
# phylDiag with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiag.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 -filterType=InFamilies -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -pThreshold=0.001 ${ibwg}identifyBreakpointsWithinGaps -overlapMax=${om} +verbose > res/`sbFileName ${S1} ${S2} ${A} ${tgm} ${gm} ${dm} ${ibwg} ${om}`"
# view the sbs from the former output file in the X-X comparison
"src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 X:1-~ X:1-~ -in:SyntenyBlocks=res/`sbFileName ${S1} ${S2} ${A} ${tgm} ${gm} ${dm} ${ibwg} ${om}` -filterType=InFamilies -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -pThreshold=0.001 ${ibwg}identifyBreakpointsWithinGaps -overlapMax=${om} +verbose -out:SyntenyBlocks=res/viewer_`sbFileName ${S1} ${S2} ${A} ${tgm} ${gm} ${dm} ${ibwg} ${om}` -out:ImageName=res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg"
)
for line in "${phylDiagCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

######################################################
#	Check integrity of phylDiagHomologyMatrixViewer.py #
######################################################

S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
gm=5
dm=CD
phylDiagHomologyMatrixViewerCommandLines=(
# phylDiagHomologyMatrixViewer -- MH
"src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 X:1-~ X:1-~ -gapMax=${gm} -distanceMetric=${dm} -out:ImageName=res/MH_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg -out:SyntenyBlocks=res/syntenyBlocksDrawerMH_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.txt"
# phylDiagHomologyMatrixViewer -- MHP
"src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 X:1-~ X:1-~ -gapMax=${gm} -distanceMetric=${dm} +mode:chromosomesRewrittenInTbs -out:ImageName=res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg -out:SyntenyBlocks=res/syntenyBlocksDrawerMHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.txt"
)
tgm=9
gm=10
dm=CD
ibwg='+'
om=10
phylDiagHomologyMatrixViewerCommandLines+=(
# phylDiagHomologyMatrixViewer with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 X:1-~ X:1-~ -tandemGapMax=${tgm} -gapMax=${gm} -distanceMetric=${dm} +mode:chromosomesRewrittenInTbs ${ibwg}identifyBreakpointsWithinGaps -overlapMax=${om} -out:ImageName=res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg -out:SyntenyBlocks=res/syntenyBlocksDrawerMHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}IbwgOm${om}.txt"
)
Title=PhylDiag
gm=5
C1=X
R1="100-210"
C2=X
R2="1-100"
phylDiagHomologyMatrixViewerCommandLines+=(
# phylDiagHomologyMatrixViewer with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${dm} -gapMax=${gm} ${C1}:${R1} ${C2}:${R2} -out:ImageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${dm}${gm}_MHP.svg -out:SyntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${dm}${gm}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001 ${ibwg}identifyBreakpointsWithinGaps -overlapMax=${om}"
#many computations known to be difficult
"src/postprocessing/listOfDifficultSyntenies.sh"
)
for line in "${phylDiagHomologyMatrixViewerCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

# Check the graph construction of the distribution of the lengths of synteny
# blocks
commandLines=(
"src/postprocessing/distributionSbLengths.py res/Homo.sapiens_Mus.musculus_Tgm10gM5Gmmi0IbwgOm10.sbs data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 -lengthUnit=Mb -minShownLength=1 -maxShownLength=81 > res/distributionOfSbsLengthsInMb.svg"
)
for line in "${commandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

# needs numpy and matplotlib


###################################
#	Check integrity of geneTeams.py #
###################################
#TODO

####################################################
#	Check integrity of geneTeamsHomologyMatrixViewer #
####################################################
Title=GeneTeams
S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
gm=4
C1=X
R1="1-~"
C2=X
R2="1-~"
geneTeamsCommandLines=(
"src/geneTeams.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 -filterType=InFamilies -gapMax=${gm} +verbose > res/geneTeams.txt"
)
for line in "${geneTeamsCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

####################################################
#	Check integrity of geneTeamsHomologyMatrixViewer #
####################################################
geneTeamsHomologyMatrixViewerCommandLines=(
"src/geneTeamsHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} -gapMax=${gm} -out:ImageName=res/${Title}_${S1}_${C1}_${S2}_${C2}_gM${gm}_MH.svg -out:GeneTeams=res/geneTeamsHomologyMatrixViewer.txt"
)
for line in "${geneTeamsHomologyMatrixViewerCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done



