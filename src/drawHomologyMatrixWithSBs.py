#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

"""
Build the homology matrix with synteny blocks using the mySvgDrawer Library
"""

import collections
import itertools
import random

import utils.mySvgDrawer as svgDrw
from utils.mySvgDrawer import Point
import utils.myLightGenomes as myLightGenomes
import utils.myGenomesDrawer as myGenomesDrawer

def test(arguments):
    scenario = arguments["scenario"]
    outFileName = arguments['out:fileName']
    if scenario==1 :
        nx=11
        ny=7
        homologySameStrand=[(1,1), (2,2), (3,3), (4,4), (5,5), (8,3), (9,4)]
        homologyOppositeStrand=[]
        diagsIndices=[[(1,1), (2,2), (3,3), (4,4), (5,5)],[(8,3), (9,4)]]
    elif scenario ==2:
        nx=9
        ny=9
        homologySameStrand=[(1,1), (2,2), (3,4), (3,5), (4,3), (5,3), (6,6), (7,7)] # not necessary to mark all for hps
        homologyOppositeStrand=[(3,3), (4,4), (4,5), (5,4), (5,5)]
        diagsIndices=[homologySameStrand + homologyOppositeStrand]
    elif scenario == 3:
        nx=6
        ny=6
        homologySameStrand=[(1,1), (2,2)]
        homologyOppositeStrand=[(3,3), (4,4)]
        diagsIndices=[homologySameStrand, [(3,3)], [(4,4)]]
    elif scenario == 4:
        nx=10
        ny=7
        homologySameStrand=[(3,3)]
        homologyOppositeStrand=[(1,5), (2,4), (4,3), (5,3), (6,3), (7, 2), (8, 1)]
        diagsIndices=[homologySameStrand + homologyOppositeStrand]
    elif scenario == 5:
        nx=9
        ny=9
        homologySameStrand=[(1,1), (2,2), (3,3), (6,6), (7,7)]
        homologyOppositeStrand=[]
        diagsIndices=[homologySameStrand, homologyOppositeStrand]
    elif scenario == 6:
        nx=13
        ny=13
        homologySameStrand=[(1,6), (2,7), (4,8), (5,9), (6,10), (10,8), (11,9)]
        homologyOppositeStrand=[(7,3),(8,2), (3,8)]
    elif scenario == 7:
        # For the paper on the method to extract diagonals
        nx=8
        ny=9
        homologySameStrand=[(1,2), (2,1), (3,3), (4,4), (4,5), (4,6), (5,4), (5,5), (5,6), (6,7)]
        homologyOppositeStrand=[(1,1),(2,2)]
        diagsIndices=[]
    elif scenario == 8:
        # For the paper on the method to extract diagonals
        nx=6
        ny=6
        homologySameStrand=[(1,1), (2,2), (3,3), (4,4)]
        homologyOppositeStrand=[]
        diagsIndices=[]
    elif scenario == 9:
        # For the paper on the method to extract diagonals
        nx=7
        ny=5
        homologySameStrand=[(2,2)]
        homologyOppositeStrand=[(1,3),(3,2),(4,2),(5,1)]
        diagsIndices=[]
    elif scenario == 10:
        # For the paper on the method to extract diagonals
        nx=5
        ny=5
        homologySameStrand=[]
        homologyOppositeStrand=[(1,3),(2,2),(3,1)]
        diagsIndices = []
    elif scenario == 11:
        # For the paper on the method to extract diagonals (merge and probabilities)
        nx=11
        ny=8
        homologyPlusAssociatedValues=[(1,(5,))]
        homologyMinusAssociatedValues=[(4,(4,)),(5,(3,)),(9,(1,))]
        homologyUnknownAssociatedValues=[(7,(6,)),(8,(2,))]
        hpSigns = collections.defaultdict(dict)
        genesStrandsC1=[]
        genesStrandsC2=[]
        for i1 in range(nx):
            genesStrandsC1.append(random.choice([-1,1]))
        for i2 in range(ny):
            genesStrandsC2.append(random.choice([-1,1]))

        for (indexesAVs,aV) in itertools.izip([homologyPlusAssociatedValues, homologyMinusAssociatedValues, homologyUnknownAssociatedValues],[+1,-1,None]):
            for (i1,i2s) in indexesAVs:
                for i2 in i2s:
                    hpSigns[i1][i2]=aV
                    # Determine possible orientations
                    strand = random.choice([-1,+1])
                    if  aV == +1:
                        genesStrandsC1[i1]=strand
                        genesStrandsC2[i2]=strand
                    elif aV == -1:
                        genesStrandsC1[i1]=strand
                        genesStrandsC2[i2]=-strand
                    elif aV == None:
                        b=random.choice([True, False])
                        genesStrandsC1[i1]=None if b else strand
                        genesStrandsC2[i2]=strand if b else None
                    else:
                        raise ValueError()
        begC1, endC1 = 0, nx-1
        begC2, endC2 = 0, ny-1
        tbWithNoHomologyInWindowC1=[0,2,3,6,10]
        tbWithNoHomologyInWindowC2=[0,7]
        homologyGroupsInWindow=[([1],[5]),([4],[4]),([5],[3]),([9],[1]),([7],[6]),([8],[2])]
        diagsIndices=[[(4,4),(5,3),(8,2),(9,1)]]
        genesRemovedDuringFilteringC1=[]
        genesRemovedDuringFilteringC2=[]
        symbolsInGenes=[[1,1,1,1,1,1,1,2,2,1,1],[1,1,2,1,1,1,2,1]]

        myGenomesDrawer.drawHomologyMatrix(((begC1, endC1),(begC2, endC2)), (genesStrandsC1, genesStrandsC2),
                           (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
                           (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2),
                           hpSigns, homologyGroupsInWindow, diagsIndices,
                           outputFileName=outFileName, maxWidth=100, maxHeight=100 ,
                           symbolsInGenes=symbolsInGenes)

    elif scenario == 12:
        genome = myLightGenomes.LightGenome()
        OG = myLightGenomes.OGene
        genome['0'] = [OG('a', +1), OG('b', +1), OG('c', -1), OG('d', +1), OG('e', -1), OG('f', +1), OG('g', +1), OG('h', +1)]
        lengthGene = 1
        width = len(genome['0']) + 2 * lengthGene
        height = 3 * lengthGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        families = myLightGenomes.Families()
        F = myLightGenomes.Family
        families.addFamily(F('A', ['a', 'b']))
        families.addFamily(F('B', ['c']))
        families.addFamily(F('C', ['d', 'e']))
        families.addFamily(F('D', ['f', 'g', 'h']))
        chromosomeItems = myGenomesDrawer.drawChromFromLightGenome(genome, '0', families, lengthGene=lengthGene)
        for item in chromosomeItems:
            # if isinstance(item, svgDrw.Gene):
            # g = item
            item.start = Point(item.start.x, item.start.y + lengthGene)
            item.end = Point(item.end.x, item.end.y + lengthGene)
            scene.add(item)
        scene.write_svg(filename=outFileName)
        print >> sys.stderr, "finished!"
    elif scenario == 13:
        families = myLightGenomes.Families()
        F = myLightGenomes.Family
        familiesDict = collections.defaultdict(list)
        allFamilies = 'ABCDEFGHIJK'
        for familyName in allFamilies:
            familiesDict[familyName] = []
        genome = myLightGenomes.LightGenome()
        OG = myLightGenomes.OGene
        print familiesDict
        for geneName in 'abcdefghijklmnopqrstuvw':
            genome['0'].append(OG(geneName, random.choice([+1, -1])))
            famName = random.choice(allFamilies)
            familiesDict[famName].append(geneName)
        for famName, dn in familiesDict.iteritems():
            families.addFamily(F(famName, dn))

        lengthGene = 1
        width = len(genome['0']) + 2 * lengthGene
        height = 3 * lengthGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        chromosomeItems = myGenomesDrawer.drawChromFromLightGenome(genome, '0', families, lengthGene=lengthGene)
        for item in chromosomeItems:
            # if isinstance(item, svgDrw.Gene):
            # g = item
            item.start = Point(item.start.x, item.start.y + lengthGene)
            item.end = Point(item.end.x, item.end.y + lengthGene)
            scene.add(item)
        scene.write_svg(filename=outFileName)
        print >> sys.stderr, "finished!"
    elif scenario == 14:
        os.chdir('/home/jlucas/Libs/MagSimus')
        genome1 = myLightGenomes.LightGenome('data/genesST.Mus.musculus.list.bz2')
        genome2 = myLightGenomes.LightGenome('data/genesST.Gallus.gallus.list.bz2')
        families = myLightGenomes.Families('data/ancGenes.Amniota.list.bz2')


        genome_Mouse = myLightGenomes.LightGenome()
        genome_Chicken = myLightGenomes.LightGenome()
        genome_Mouse['5'] = genome1['5'][597-1:627-1]
        genome_Chicken['4'] = genome2['4'][582-1:611-1]

        newFamilies = myLightGenomes.Families()
        old2newFamily = {}
        allFamilies = []
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789':
            allFamilies.append(c)

        familiesInBothGenomes = genome_Mouse.getOwnedFamilyNames(families) & genome_Chicken.getOwnedFamilyNames(families)
        for gene in genome_Mouse['5'] + genome_Chicken['4']:
            oldFamily = families.getFamilyByName(gene.n, default=None)
            if oldFamily:
                if oldFamily.fn in familiesInBothGenomes:
                    if oldFamily.fn not in old2newFamily:
                        newFamily = allFamilies.pop(0)
                        newFamilies.addFamily(myLightGenomes.Family(newFamily, oldFamily.dns))
                        old2newFamily[oldFamily.fn] = newFamily
        families = newFamilies

        lengthGene = 1
        width = max(len(genome_Mouse['5']), len(genome_Chicken['4'])) + 2 * lengthGene
        height = 4 * lengthGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        chromosomeMouseItems = myGenomesDrawer.drawChromFromLightGenome(genome_Mouse, '5', families, lengthGene=lengthGene)
        chromosomeChickenItems = myGenomesDrawer.drawChromFromLightGenome(genome_Chicken, '4', families, lengthGene=lengthGene)
        for item in chromosomeMouseItems:
            item.start = Point(item.start.x, item.start.y + lengthGene)
            item.end = Point(item.end.x, item.end.y + lengthGene)
            scene.add(item)
        for item in chromosomeChickenItems:
            item.start = Point(item.start.x, item.start.y + 2*lengthGene)
            item.end = Point(item.end.x, item.end.y + 2*lengthGene)
            scene.add(item)
        scene.write_svg(filename=outFileName)

    elif scenario == 15:
        # presentation of the different events
        OG = myLightGenomes.OGene
        genome = myLightGenomes.LightGenome()
        genome['1'] = [OG('A', +1), OG('B', -1), OG('C', +1)]
        genome['2'] = [OG('D', -1), OG('E', +1), OG('F', -1)]

        sizeGene = 1
        listOfChromosomes = []
        for c in genome:
            listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(genome, c, lengthGene=sizeGene, 
                                                                                     homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        import libs.myEvents as mE

        # Gene events
        # tandem duplication
        newGenome = mE.performInsertNewGene(genome, ('B.b', '1', 2, -1), keepOriginalGenome=True)
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        # gene loss
        newGenome = mE.performGeneLoss(genome, ('1', 1), keepOriginalGenome=True)
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        # de novo gene birth
        newGenome = mE.performInsertNewGene(genome, ('G', '1', 2, -1), keepOriginalGenome=True)
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))

        # Chromosomal rearrangements
        # fission
        newGenome = mE.performFission(genome, ('1', 1), keepOriginalGenome=True)
        assert '0' in newGenome.keys()
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '0', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        # fusion
        newGenome = mE.performFusion(genome, (('1', +1), ('2', +1)), keepOriginalGenome=True)
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        # inversion
        newGenome = mE.performInversion(genome, ('1', 1, 3), keepOriginalGenome=True)
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        # reciprocal translocation
        newGenome = mE.performReciprocalTranslocation(genome, (('1', 1, +1), ('2', 1, +1)), keepOriginalGenome=True)
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '1', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(newGenome, '2', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))

        genome['3'] = [OG('A', +1), OG('B', -1), OG('C', +1), OG('D', +1), OG('E', -1), OG('F', +1), OG('G', +1), OG('H', -1), OG('I', +1)]
        listOfChromosomes.append(myGenomesDrawer.drawChromFromLightGenome(genome, '3', lengthGene=sizeGene, homologsColorsGenerator=myGenomesDrawer.levelIdxGenerator(farIdxs=5)))

        width = (2 + 8) * sizeGene
        height = (2 + len(listOfChromosomes)) * sizeGene
        scene = svgDrw.Scene(name='chromosome', width=width, height=height)
        for i, chromosome in enumerate(listOfChromosomes):
            svgDrw.tanslateItems(chromosome, 0, sizeGene + i * sizeGene)
            for item in chromosome:
                scene.add(item)
        scene.write_svg(filename=outFileName)



if __name__ == '__main__':
    #arguments = myTools.checkArgs([("scenario",int)],[("out:FileName",str,"image.svg")],__doc__)
    import os
    import sys
    print sys.stderr, sys.argv
    os.chdir('/home/jlucas/Libs/MagSimus')
    arguments = {}
    arguments['scenario'] = int(sys.argv[1])
    arguments['out:fileName'] = 'toto.svg'
    test(arguments)
