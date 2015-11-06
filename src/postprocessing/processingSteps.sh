#!/bin/bash
#Data78
set -e
green='\e[0;32m'
NC='\e[0m' # No Color

# scale factor for rectangles
sffr=5
pt=1
dm=CD

# Human-Mouse comparison, the tandem duplication solved thanks to the
# truncation and re-merge
S1=Homo.sapiens
S2=Mus.musculus
A=Euarchontoglires
C1=X
beg1=1
end1='~'
C2=X
beg2=1
end2='~'
R1="${beg1}-${end1}"
R2="${beg2}-${end2}"

mkdir res/processingSteps &

cmd=()
# MH no filtering
ritb='-'
tgm=0
gm=0
ibwg='-'
no='-'
om=0
filter='None'
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:ImageName=res/processingSteps/MH_0.svg -out:SyntenyBlocks=res/processingSteps/MH_0.txt -scaleFactorRectangles=${sffr}")

# tandemGapMax filtering InFamilies
filter='InBothGenomes'
for tgm in 0 10
do
	cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:ImageName=res/processingSteps/MH_InFamilies_tgm${tgm}.svg -out:SyntenyBlocks=res/processingSteps/MH_InFamilies_tgm${tgm}.txt -scaleFactorRectangles=${sffr}")
done

# tandemGapMax filtering InBothGenomes
#filter='InBothGenomes'
#tgm=10
#cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:ImageName=res/processingSteps/MH_InBothGenomes_tgm${tgm}.svg -out:SyntenyBlocks=res/processingSteps/MH_InBothGenomes_tgm${tgm}.txt -scaleFactorRectangles=${sffr}")

# gapMaxs
for gm in 1 2 3 4 5 10
do
	cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:ImageName=res/processingSteps/MH_InBothGenomes_tgm${tgm}_gm${gm}.svg -out:SyntenyBlocks=res/processingSteps/MH_InBothGenomes_tgm${tgm}_gm${gm}.txt -scaleFactorRectangles=${sffr}")
done

# identify breakpoints within gaps
ibwg='+'
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:ImageName=res/processingSteps/MH_InBothGenomes_tgm${tgm}_gm${gm}_ibwg.svg -out:SyntenyBlocks=res/processingSteps/MH_InBothGenomes_tgm${tgm}_gm${gm}_ibwg.txt -scaleFactorRectangles=${sffr}")

# no overlap
no='+'
# overlap maxs
for om in 0 5 10
do
	cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:ImageName=res/processingSteps/MH_InBothGenomes_tgm${tgm}_gm${gm}_ibwg_om${om}.svg -out:SyntenyBlocks=res/processingSteps/MH_InBothGenomes_tgm${tgm}_gm${gm}_ibwg_om${om}.txt -scaleFactorRectangles=${sffr}")
done

for line in "${cmd[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

