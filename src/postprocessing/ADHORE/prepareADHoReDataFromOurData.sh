#!/bin/bash

# ligne de commande classique
# ./prepareADHoReDataFromMagSimus.sh ../../PhylDiag/data/genesSTE.Homo.sapiens.list.bz2 ../../PhylDiag/data/genesSTE.Mus.musculus.list.bz2 ../../PhylDiag/data/ancGenes.Euarchontoglires.list.bz2 5 Hs_Mm

# For condor over the 100 simulations
# ../../MagSimus/res/simu1/
genome1=$1
genome2=$2
family=$3
maxGap=$4
outDir=$5_${maxGap}
mkdir -p ./${outDir}/G1/
./genes2ADHoReData.py ${genome1} ${family} -out:Chromosomes=${outDir}/G1/Genome.G1.Chr%s.list > ${outDir}/families.csv 2> ${outDir}/G1/logConversion1

mkdir -p ./${outDir}/G2/
./genes2ADHoReData.py ${genome2} ${family} -out:Chromosomes=${outDir}/G2/Genome.G2.Chr%s.list >> ${outDir}/families.csv 2> ${outDir}/G2/logConversion2
sort -k2 ${outDir}/families.csv > ${outDir}/families.csv_
mv ${outDir}/families.csv_ ${outDir}/families.csv
# copy the configuration file
cp ./dataset_G1_G2.ini ${outDir}/dataset_G1_G2_DPD${maxGap}.ini
# Configure maxGap parameter
sed -i "/output_path/s/=.*/=OutPutADHoRe_DPD${maxGap}/" ${outDir}/dataset_G1_G2_DPD${maxGap}.ini
sed -i "/gap_size/s/=.*/=${maxGap}/" ${outDir}/dataset_G1_G2_DPD${maxGap}.ini
sed -i "/cluster_gap/s/=.*/=${maxGap}/" ${outDir}/dataset_G1_G2_DPD${maxGap}.ini
