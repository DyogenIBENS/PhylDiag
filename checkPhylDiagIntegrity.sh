#!/bin/bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

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
	# identify micro-rearrangements within gaps
	imrwg=$7
	# truncation max
	tm=$8
	if [ "${imrwg}" == "+" ]; then
		echo "${s1}_${s2}_f${a}_Tgm${tgm}Gm${gm}${dm}ImrwgTm${om}.sbs"
	else
		echo "${s1}_${s2}_f${a}_Tgm${tgm}Gm${gm}${dm}Tm${om}.sbs"
	fi
}

#############################################
#	Check integrity of phylDiag.py					  #
#############################################
data="data/"
S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
phylDiagEasyCommandLines=(
"src/phylDiag.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 > res/${S1}_${S2}_f${A}.cs"
)
S1=Homo.sapiens
S2=Gallus.gallus
A=Amniota
phylDiagEasyCommandLines+=(
"src/phylDiag.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 > res/${S1}_${S2}_f${A}.cs"
)
f='InBothGenomes'
tgm=9
gm=5
dm=CD
#imrwg
tm=10
mmg=1
#verbose
phylDiagEasyCommandLines+=(
"src/phylDiag.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 -f $f -t ${tm} -d ${dm} -g ${gm} --imr --imcs --mmg ${mmg} --truncation --truncationMax ${tm} --verbose > res/${S1}_${S2}_f${A}.cs"
)
f='InBothGenomes'
tgm=3
gm=20
dm=ED
#imrwg
tm=2
mmg=3
phylDiagEasyCommandLines+=(
"src/phylDiag.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 -f $f -t ${tm} -d ${dm} -g ${gm} --imr --imcs --mmg ${mmg} --truncation --truncationMax ${tm} --verbose > res/${S1}_${S2}_f${A}.cs"
)

for line in "${phylDiagEasyCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done



#############################################
#	Check integrity of phylDiagCA.py					  #
#############################################

S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
phylDiagCommandLines=(
# phylDiagCA with default options
"src/phylDiagCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 > res/${S1}_${S2}_f${A}.sbs"
)
gm=5
dm=CD
phylDiagCommandLines+=(
# phylDiag with some options
"src/phylDiagCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 -filterType=InFamilies -distanceMetric=${dm} -gapMax=${gm} > res/${S1}_${S2}_f${A}_Gm${gm}${dm}.sbs"
"src/phylDiagCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 -filterType=InBothGenomes -distanceMetric=${dm} -gapMax=${gm} > res/${S1}_${S2}_fIBS${A}_Gm${gm}${dm}.sbs"
)
tgm=9
gm=5
dm=CD
imrwg='+'
tm=10
phylDiagCommandLines+=(
# phylDiag with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiagCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 -filterType=InFamilies -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -pThreshold=0.001 ${imrwg}identifyMicroRearr -truncationMax=${tm} +verbose > res/`sbFileName ${S1} ${S2} ${A} ${tgm} ${gm} ${dm} ${imrwg} ${tm}`"
# view the sbs from the former output file in the X-X comparison
"src/phylDiagViewerCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 X:1-~ X:1-~ -in:syntenyBlocks=res/`sbFileName ${S1} ${S2} ${A} ${tgm} ${gm} ${dm} ${imrwg} ${tm}` -filterType=InFamilies -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -pThreshold=0.001 ${imrwg}identifyMicroRearr -truncationMax=${tm} +verbose -out:syntenyBlocks=res/viewer_`sbFileName ${S1} ${S2} ${A} ${tgm} ${gm} ${dm} ${imrwg} ${tm}` -out:imageFileName=res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg"
)
for line in "${phylDiagCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

##################################################
#	Check integrity of phylDiagViewer.py         #
##################################################

S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
imrwg='no-' # means no
gm=5
dm=CD
phylDiagViewerCommandLines=(
# phylDiagViewer -- MH
"src/phylDiagViewer.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 --ROI1=X:1-~ --ROI2=X:1-~ --gapMax=${gm} --distanceMetric=${dm} --${imrwg}imr res/MH_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg --outSbs=res/syntenyBlocksDrawerMH_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.txt"
# phylDiagViewer -- MHP
"src/phylDiagViewer.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 --ROI1=X:1-~ --ROI2=X:1-~ --gapMax=${gm} --distanceMetric=${dm} --chrsInTbs --${imrwg}imr res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg --outSbs=res/syntenyBlocksDrawerMHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.txt"
)
tgm=9
gm=10
dm=CD
imrwg='' # means yes
tm=10
phylDiagViewerCommandLines+=(
# phylDiagViewer with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiagViewer.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 --ROI1=X:1-~ --ROI2=X:1-~ --tandemGapMax=${tgm} --gapMax=${gm} --distanceMetric=${dm} --chrsInTbs --${imrwg}imr --truncationMax=${tm} res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg --outSbs=res/syntenyBlocksDrawerMHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}ImrwgTm${tm}.txt"
)
Title=PhylDiag
gm=5
C1=X
R1="100-210"
C2=X
R2="1-100"
phylDiagViewerCommandLines+=(
# phylDiagViewer with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiagViewer.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 --chrsInTbs --distanceMetric=${dm} --gapMax=${gm} --ROI1=${C1}:${R1} --ROI2=${C2}:${R2} res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${dm}${gm}_MHP.svg --outSbs=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${dm}${gm}_syntenyBlocksDrawerMHP.txt --verbose --${imrwg}imr --truncationMax=${tm} --outSbs=res/syntenyBlocksDrawer.txt"
#many computations known to be difficult
"src/postprocessing/listOfDifficultSyntenies.sh"
)
for line in "${phylDiagViewerCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

# TODO
# needs numpy and matplotlib
# Check the graph construction of the distribution of the lengths of synteny
# blocks
#commandLines=(
#"src/postprocessing/distributionSbLengths.py res/Homo.sapiens_Mus.musculus_Tgm10gM5Gmmi0ImrwgTm10.sbs data/Homo.sapiens.genome.bz2 data/Mus.musculus.genome.bz2 -lengthUnit=Mb -minSbLength=1 -maxSbLength=81 > res/distributionOfSbsLengthsInMb.svg"
#)
#for line in "${commandLines[@]}"
#	do
#		echo -e "${green}${line}${NC}"
#		eval ${line}
#done


######################################################
#	Check integrity of phylDiagViewerCA.py #
######################################################

S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
gm=5
dm=CD
phylDiagViewerCACommandLines=(
# phylDiagViewerCA -- MH
"src/phylDiagViewerCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 X:1-~ X:1-~ -gapMax=${gm} -distanceMetric=${dm} -out:imageFileName=res/MH_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg -out:syntenyBlocks=res/syntenyBlocksDrawerMH_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.txt"
# phylDiagViewerCA -- MHP
"src/phylDiagViewerCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 X:1-~ X:1-~ -gapMax=${gm} -distanceMetric=${dm} +mode:chromosomesRewrittenInTbs -out:imageFileName=res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg -out:syntenyBlocks=res/syntenyBlocksDrawerMHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.txt"
)
tgm=9
gm=10
dm=CD
imrwg='+'
tm=10
phylDiagViewerCACommandLines+=(
# phylDiagViewerCA with all options -- options that yield the best synteny blocks for the human mouse comparison
"src/phylDiagViewerCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 X:1-~ X:1-~ -tandemGapMax=${tgm} -gapMax=${gm} -distanceMetric=${dm} +mode:chromosomesRewrittenInTbs ${imrwg}identifyMicroRearr -truncationMax=${tm} -out:imageFileName=res/MHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}.svg -out:syntenyBlocks=res/syntenyBlocksDrawerMHP_${S1}_X_${S2}_X_f${A}_Gm${gm}${dm}ImrwgTm${om}.txt"
)
Title=PhylDiag
gm=5
C1=X
R1="100-210"
C2=X
R2="1-100"
phylDiagViewerCACommandLines+=(
# phylDiagViewerCA with all options -- options that yield the best synteny blocks for the human mouse comparison
# FIXME -out:syntenyBlocks=res/syntenyBlocksDrawer
"src/phylDiagViewerCA.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${dm} -gapMax=${gm} ${C1}:${R1} ${C2}:${R2} -out:imageFileName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${dm}${gm}_MHP.svg -out:syntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${dm}${gm}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001 ${imrwg}identifyMicroRearr -truncationMax=${tm} -out:syntenyBlocks=res/syntenyBlocksDrawer.txt"
#many computations known to be difficult
"src/postprocessing/listOfDifficultSyntenies.sh"
)
for line in "${phylDiagViewerCACommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

# TODO
# needs numpy and matplotlib
# Check the graph construction of the distribution of the lengths of synteny
# blocks
#commandLines=(
#"src/postprocessing/distributionSbLengths.py res/Homo.sapiens_Mus.musculus_Tgm10gM5Gmmi0ImrwgTm10.sbs data/Homo.sapiens.genome.bz2 data/Mus.musculus.genome.bz2 -lengthUnit=Mb -minSbLength=1 -maxSbLength=81 > res/distributionOfSbsLengthsInMb.svg"
#)
#for line in "${commandLines[@]}"
#	do
#		echo -e "${green}${line}${NC}"
#		eval ${line}
#done

#####################################################
#	Check integrity of wholeGenomeHMCA.py           #
#####################################################
commandLines=(
"src/wholeGenomeHMCA.py data/Homo.sapiens.genome.bz2 data/Mus.musculus.genome.bz2 data/Euarchontoglires.families.bz2 -tandemGapMax=5 -gapMax=5 -out:imageFileName=res/WMH_Hs_Mm.svg +withSbs -scaleFactorRectangles=40"
)
for line in "${commandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

#########################################
#	Check integrity of geneTeams.py #
#########################################
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
"src/geneTeams.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 -filterType=InFamilies -gapMax=${gm} +verbose > res/geneTeams.txt"
)
for line in "${geneTeamsCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done


##########################################################
#	Check integrity of the benchmark scripts         #
##########################################################
data="data/benchmark"
benchmarkCmdLines=(
"src/benchmark/comparePhylDiagAndADHORESbsToSimulatedSbs.py ${data}/speciesTree.phylTree Mus.musculus Gallus.gallus Amniota -pSimGenomes=${data}/%s.genome.bz2 -pAncGenes=${data}/%s.families.bz2 -pSimulatedSbs=${data}/%s.%s.sbs.bz2  -preComputePairwiseSbs -oriented")
benchmarkCmdLines+=(
"src/benchmark/compareCyntenatorSbsToSimulatedSbs.py ${data}/speciesTree.phylTree Mus.musculus Gallus.gallus Amniota -pSimGenomes=${data}/%s.genome.bz2 -pAncGenes=${data}/%s.families.bz2 -pSimulatedSbs=%s.%s.sbs.bz2 -preComputePairwiseSbs -oriented")
for line in "${benchmarkCmdLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

# TODO
##########################################################
#	Check integrity of geneTeamsHomologyMatrixViewer #
##########################################################
#geneTeamsHomologyMatrixViewerCommandLines=(
#"src/geneTeamsHomologyMatrixViewer.py data/${S1}.genome.bz2 data/${S2}.genome.bz2 data/${A}.families.bz2 ${C1}:${R1} ${C2}:${R2} -gapMax=${gm} -out:imageName=res/${Title}_${S1}_${C1}_${S2}_${C2}_gM${gm}_MH.svg -out:geneTeams=res/geneTeamsHomologyMatrixViewer.txt"
#)
#for line in "${geneTeamsHomologyMatrixViewerCommandLines[@]}"
#	do
#		echo -e "${green}${line}${NC}"
#		eval ${line}
#done



