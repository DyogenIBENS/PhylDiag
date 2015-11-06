#!/bin/bash
#Data78
set -e
green='\e[0;32m'
NC='\e[0m' # No Color

# scale factor for rectangles
sffr=5

tgm=9
gm=10
dm=CD
ibwg='+'
om=10
# p-value threshold
pt=1

# Human-Mouse comparison, the tandem duplication solved thanks to the
# truncation and re-merge
S1=Mus.musculus
S2=Gallus.gallus
A=Amniota
C1=5
beg1=597
end1=626
C2=4
beg2=582
end2=610
R1="${beg1}-${end1}"
R2="${beg2}-${end2}"

mkdir res/processingSteps_Mm_Gg &

cmd=()
# MH no filtering
ritb='-'
tgm=0
gm=0
ibwg='-'
no='-'
om=1000
filter='None'
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -filterType=${filter} -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MH_0.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MH_0.txt -scaleFactorRectangles=${sffr}")

# tandemGapMax filtering InFamilies
filter='InFamilies'
ritb='+'
gm=0
ibwg='-'
no='-'
om=0
for tgm in 0 10
do
	cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InFamilies_tgm${tgm}.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InFamilies_tgm${tgm}.txt -scaleFactorRectangles=${sffr} +convertGenicToTbCoordinates")
done

# tandemGapMax filtering InBothGenomes
filter='InBothGenomes'
tgm=10
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}.txt -scaleFactorRectangles=${sffr} +convertGenicToTbCoordinates")

# gapMaxs
for gm in 0 5 10
do
	cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}.txt -scaleFactorRectangles=${sffr} +convertGenicToTbCoordinates")
done

# identify breakpoints within gaps
ibwg='+'
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg.txt -scaleFactorRectangles=${sffr} +convertGenicToTbCoordinates")

# no overlap
no='+'
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg.txt -scaleFactorRectangles=${sffr} +convertGenicToTbCoordinates")

# overlap maxs
for om in 0 5 10
do
	cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg_om${om}.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg_om${om}.txt -scaleFactorRectangles=${sffr} +convertGenicToTbCoordinates")
done
# last with MH
ritb='-'
cmd+=("src/phylDiagHomologyMatrixViewer.py data/genesST.${S1}.list.bz2 data/genesST.${S2}.list.bz2 data/ancGenes.${A}.list.bz2 ${C1}:${R1} ${C2}:${R2} ${ritb}mode:chromosomesRewrittenInTbs -tandemGapMax=${tgm} -distanceMetric=${dm} -gapMax=${gm} -verbose -pThreshold=${pt} ${ibwg}identifyMicroRearrangements ${no}nonOverlappingSbs -truncationMax=${om} -out:imageName=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg_om${om}.svg -out:syntenyBlocks=res/processingSteps_Mm_Gg/MHP_InBothGenomes_tgm${tgm}_gm${gm}_ibwg_om${om}.txt -scaleFactorRectangles=${sffr}")

for line in "${cmd[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done

