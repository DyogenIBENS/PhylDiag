#!/bin/bash
#Data78
set -e
green='\e[0;32m'
NC='\e[0m' # No Color

# Human-Mouse comparison
args=(
#S1	S2	LCA	C1 beg1 end1 C2 beg2 end2
"Homo.sapiens Mus.musculus Euarchontoglires 1 1 ~ 4 1 ~"
"Homo.sapiens Mus.musculus Euarchontoglires 6 1 ~ 17 1 ~"
"Homo.sapiens Mus.musculus Euarchontoglires 17 1 ~ 11 1 ~"
"Homo.sapiens Mus.musculus Euarchontoglires 19 1 ~ 7 1 ~"
"Homo.sapiens Mus.musculus Euarchontoglires 20 1 ~ 2 1 ~"
"Homo.sapiens Mus.musculus Euarchontoglires X 1 ~ X 1 ~"
)
# Human-Chicken comparison
args+=(
#S1	S2	LCA	C1 beg1 end1 C2 beg2 end2
"Homo.sapiens Gallus.gallus Amniota 6 1 ~ 3 1 ~"
"Homo.sapiens Gallus.gallus Amniota X 1 ~ 4 1 ~"
"Homo.sapiens Gallus.gallus Amniota 4 1 ~ 4 1 ~"
"Homo.sapiens Gallus.gallus Amniota 10 1 ~ 6 1 ~"
"Homo.sapiens Gallus.gallus Amniota 16 1 ~ 11 1 ~"
)
# Mouse-Chicken comparison
args=(
#S1	S2	LCA	C1 beg1 end1 C2 beg2 end2
"Mus.musculus Gallus.gallus Amniota X 1 ~ 4 1 ~"
"Mus.musculus Gallus.gallus Amniota 5 1 ~ 4 1 ~"
)

Title=MHP
C1=X
R1="100-250"
C2=X
R2="1-100"
tgm=9
gm=10
dm=CD
ibwg='+'
om=10

commandLines=()
for line in "${args[@]}"
do

	S1=$(echo ${line}|cut -d" " -f1)
	S2=$(echo ${line}|cut -d" " -f2)
	A=$(echo ${line}|cut -d" " -f3)
	C1=$(echo $line|cut -d" " -f4)
	beg1=$(echo $line|cut -d" " -f5)
	end1=$(echo $line|cut -d" " -f6)
	C2=$(echo $line|cut -d" " -f7)
	beg2=$(echo $line|cut -d" " -f8)
	end2=$(echo $line|cut -d" " -f9)
	R1="${beg1}-${end1}"
	R2="${beg2}-${end2}"

	commandLines+=(
	"src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} +mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=0.001 ${ibwg}identifyBreakpointsWithinGaps +nonOverlappingSbs -overlapMax=${om} -out:ImageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}.svg -out:SyntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}_syntenyBlocksDrawer.txt"
	)
done

for line in "${commandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done
