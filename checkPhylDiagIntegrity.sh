#!/bin/bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent commands which fail will cause the shell script to exit
# immediately
red='\e[0;31m'
NC='\e[0m' # No Color

# TODO for each command line, avoid duplicates
#commandLine = "cmd"
##echo -e "${cmd}${NC}"
#${cmd}

#echo -e "${red}src/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.protTree${NC}"
#src/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.protTree
#echo -e "${red}src/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree${NC}"
#src/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree
#echo -e "${red}src/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.protTree${NC}"
#src/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.protTree
#
#echo -e "${red}src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -filterType=InCommonAncestor > res/syntenyBlocks.txt${NC}"
#src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -filterType=InCommonAncestor > res/syntenyBlocks.txt
#echo -e "${red}src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -gapMax=11 -distanceMetric=MD -filterType=InCommonAncestor > res/syntenyBlocks_MD11.txt${NC}"
#src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -gapMax=11 -distanceMetric=MD -filterType=InCommonAncestor > res/syntenyBlocks_MD11.txt
#echo -e "${red}src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -distanceMetric=DPD -gapMax=5 -pThreshold=0.001 +verbose -filterType=InCommonAncestor > res/syntenyBlocks_DPD11.txt 2> res/logErr_DPD11.txt${NC}"
#src/phylDiag.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -distanceMetric=DPD -gapMax=5 -pThreshold=0.001 +verbose -filterType=InCommonAncestor > res/syntenyBlocks_DPD11.txt
#
#echo -e "${red}src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -gapMax=5 -distanceMetric=MD -out:ImageName=res/MH_MD5.svg -out:SyntenyBlocks=res/syntenyBlocksDrawerMH_MD5.txt${NC}"
#src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -gapMax=5 -distanceMetric=MD -out:ImageName=res/MH_MD5.svg -out:SyntenyBlocks=res/syntenyBlocksDrawerMH_MD5.txt
#echo -e "${red}src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ +mode:chromosomesRewrittenInTbs -distanceMetric=MD -gapMax=5 -out:ImageName=./res/MHP_MD5.svg -out:SyntenyBlocks=./res/syntenyBlocksDrawerMHP_MD5.txt${NC}"
#src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ +mode:chromosomesRewrittenInTbs -distanceMetric=MD -gapMax=5 -out:ImageName=./res/MHP_MD5.svg -out:SyntenyBlocks=./res/syntenyBlocksDrawerMHP_MD5.txt
#
#Title=PhylDiag
#S1=Homo.sapiens
#S2=Mus.musculus
#C1=X
#R1="100-250"
#C2=X
#R2="1-100"
#DM="DPD"
#D=10
#echo -e "${red}Title=PhylDiag && S1=Homo.sapiens && S2=Mus.musculus && C1=X && R1="100-250" && C2=X && R2="1-100" && DM="DPD" && D=10 && src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${DM} -gapMax=${D} $C1:$R1 $C2:$R2 -out:ImageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_MHP.svg -out:SyntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001${NC}"
#src/phylDiagHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 +mode:chromosomesRewrittenInTbs -distanceMetric=${DM} -gapMax=${D} ${C1}:${R1} ${C2}:${R2} -out:ImageName=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_MHP.svg -out:SyntenyBlocks=res/${Title}_${S1}_${C1}.${R1}_${S2}_${C2}.${R2}_${DM}${D}_syntenyBlocksDrawerMHP.txt -verbose -pThreshold=0.001

# check geneTeams integrity
Title=GeneTeams
echo -e "${red}src/geneTeamsHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -gapMax=3 -out:ImageName=res/${NC}"
src/geneTeamsHomologyMatrixViewer.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 X:1-~ X:1-~ -gapMax=5 -out:ImageName=res/${Title}_${S1}_${C1}_${S2}_${C2}_${D}_MH.svg -out:GeneTeams=res/geneTeams.txt
echo -e "${red}src/geneTeams.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -gapMax=5 +verbose -filterType=InCommonAncestor > res/geneTeams.txt${NC}"
src/geneTeams.py data/genesST.Homo.sapiens.list.bz2 data/genesST.Mus.musculus.list.bz2 data/ancGenes.Euarchontoglires.list.bz2 -gapMax=5 +verbose -filterType=InCommonAncestor > res/geneTeams.txt
