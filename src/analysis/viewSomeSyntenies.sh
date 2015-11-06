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
# Human-Dog comparison
#args+=(
##S1	S2	LCA	C1 beg1 end1 C2 beg2 end2
#"Homo.sapiens Canis.lupus.familiaris Boreoeutheria 6 1 ~ 3 1 ~"
#"Homo.sapiens Canis.lupus.familiaris Boreoeutheria X 1 ~ 4 1 ~"
#"Homo.sapiens Canis.lupus.familiaris Boreoeutheria 4 1 ~ 4 1 ~"
#"Homo.sapiens Canis.lupus.familiaris Boreoeutheria 10 1 ~ 6 1 ~"
#"Homo.sapiens Canis.lupus.familiaris Boreoeutheria 16 1 ~ 11 1 ~"
#)
## Mouse-Cow comparison
#args+=(
##S1	S2	LCA	C1 beg1 end1 C2 beg2 end2
#"Mus.musculus Bos.taurus Boreoeutheria X 1 ~ 4 1 ~"
#"Mus.musculus Bos.taurus Boreoeutheria 5 1 ~ 4 1 ~"
#)

#Title=MHP
C1=X
#R1="100-250"
C2=X
#R2="1-100"
f=InBothGenomes
tgm=7
gm=5
dm=CD
ibwg='+'
gmmi=2
om=None
# p-value threshold
pt=None

commandLines=()

# 1: compute sbs
speciesCombi=(
"Homo.sapiens Mus.musculus Euarchontoglires"
#"Homo.sapiens Canis.lupus.familiaris Boreoeutheria"
#"Mus.musculus Bos.taurus Boreoeutheria"
)
for line in "${speciesCombi[@]}"
do
S1=$(echo ${line}|cut -d" " -f1)
S2=$(echo ${line}|cut -d" " -f2)
A=$(echo ${line}|cut -d" " -f3)
commandLines+=(
"python -O ../phylDiag.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 -filterType=${f} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -pThreshold=${pt} -gapMaxMicroInv=${gmmi} ${ibwg}identifyMicroRearrangements -truncationMax=${om} -verbose > res/${S1}_${S2}_Tgm${tgm}gM${gm}Gmmi${gmmi}IbwgOm${om}.sbs 2> >(tee logErr_syntenyBlocksDrawer.txt >&2)"
)
done

#2: show sbs for each ROI
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
	# MH
	# recompute the sbs
	#"python -O src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} -filterType=${f} -mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} -gapMaxMicroInv=${gmmi} ${ibwg}identifyMicroRearrangements -truncationMax=${om} -out:ImageName=res/MH_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}.svg -out:SyntenyBlocks=res/MH_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}_syntenyBlocksDrawer.txt"
	# reuse previously computed sbs
	"python -O ../phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} -filterType=${f} -mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} -gapMaxMicroInv=${gmmi} ${ibwg}identifyMicroRearrangements -truncationMax=${om} -out:ImageName=res/MH_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}.svg -out:SyntenyBlocks=res/MH_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}_syntenyBlocksDrawer.txt -in:SyntenyBlocks=res/${S1}_${S2}_Tgm${tgm}gM${gm}Gmmi${gmmi}IbwgOm${om}.sbs"

	# MHP
	# recompute the sbs
	#"python -O src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} -filterType=${f} +mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} -gapMaxMicroInv=${gmmi} ${ibwg}identifyMicroRearrangements -truncationMax=${om} -out:ImageName=res/MHP_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}.svg -out:SyntenyBlocks=res/MHP_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}_syntenyBlocksDrawer.txt"
	# reuse previously computedsbs
	"python -O ../phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} -filterType=${f} +mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} -gapMaxMicroInv=${gmmi} ${ibwg}identifyMicroRearrangements -truncationMax=${om} -out:ImageName=res/MHP_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}.svg -out:SyntenyBlocks=res/MHP_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_Tgm${tgm}gM${gm}IbwgOm${om}_syntenyBlocksDrawer.txt -in:SyntenyBlocks=res/${S1}_${S2}_Tgm${tgm}gM${gm}Gmmi${gmmi}IbwgOm${om}.sbs"

	)
done

for line in "${commandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done
