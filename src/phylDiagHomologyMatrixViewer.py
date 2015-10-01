#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import utils.myLightGenomes as myLightGenomes
import utils.myTools as myTools
import utils.myGenomesDrawer as myGenomesDrawer

# PhylDiag core algorithm
import utils.myDiags as myDiags

filterType = list(myDiags.FilterType._keys)

__doc__ = """
        Show the homology matrix with coloured synteny blocks (also called diagonals).
        - On the x-axis are the genes of the 1st genome in the desired window
        - On the y-axis are the genes of the 2nd genome in the desired window
        - Each coloured rectangle in the matrix represents a homology
        - '+' indicates homology that have horizontal gene and vertical gene in the same direction on both chromosomes
          '-' indicates homology that have horizontal gene and vertical gene in opposite directions on each chromosomes
        """

arguments = myTools.checkArgs(
    [("genome1", file),
     ("genome2", file),
     ("families", file),
     ("chr1:deb1-fin1", str),
     ("chr2:deb2-fin2", str)],
    myDiags.defaultArgsPhylDiag +\
    [('removeUnofficialChromosomes', bool, True),
     ('withSbs', bool, True),
     ("in:SyntenyBlocks", str, 'None'),
     ("out:SyntenyBlocks", str, "./syntenyBlocksDrawer.txt"),
     ("mode:chromosomesRewrittenInTbs", bool, False),
     ('convertGenicToTbCoordinates', bool, False),
     ('drawAllInformations', bool, False),
     ("scaleFactorRectangles", float, 2.0),
     ("out:ImageFileName", str, "./homologyMatrix.svg"),
     ("considerAllPairComps", bool, True),
     ('switchOnDirectView', bool, False),
     ('verbose', bool, True)],
    __doc__)

kwargs = myDiags.defaultKwargsPhylDiag(arguments=arguments)
kwargs['verbose'] = True

genome1 = myLightGenomes.LightGenome(arguments['genome1'], withDict=True)
genome2 = myLightGenomes.LightGenome(arguments['genome2'], withDict=True)
families = myLightGenomes.Families(arguments['families'])

if arguments['removeUnofficialChromosomes']:
    genome1.removeUnofficialChromosomes()
    genome2.removeUnofficialChromosomes()

sbsInPairComp = None
if arguments['withSbs']:
    if arguments['in:SyntenyBlocks'] != 'None':
        # load precomputed sbs if any
        sbsInPairComp = myDiags.parseSbsFile(arguments['in:SyntenyBlocks'], genome1=genome1, genome2=genome2)
    else:
        sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, families,
                                                            **kwargs)

myGenomesDrawer.homologyMatrixViewer(genome1, genome2, families, arguments['chr1:deb1-fin1'], arguments['chr2:deb2-fin2'],
                                     convertGenicToTbCoordinates=arguments['convertGenicToTbCoordinates'],
                                     filterType=kwargs['filterType'],
                                     minChromLength=kwargs['minChromLength'],
                                     tandemGapMax=kwargs['tandemGapMax'],
                                     distanceMetric=kwargs['distanceMetric'],
                                     gapMax=kwargs['gapMax'],
                                     distinguishMonoGenicDiags=kwargs['distinguishMonoGenicDiags'],
                                     pThreshold=kwargs['pThreshold'],
                                     gapMaxMicroInv=kwargs['gapMaxMicroInv'],
                                     identifyMonoGenicInversion=kwargs["identifyMonoGenicInversion"],
                                     identifyBreakpointsWithinGaps=kwargs['identifyBreakpointsWithinGaps'],
                                     overlapMax=kwargs['overlapMax'],
                                     sameStrand=kwargs['sameStrand'],
                                     validateImpossToCalc_mThreshold=kwargs['validateImpossToCalc_mThreshold'],
                                     nbHpsRecommendedGap=kwargs['nbHpsRecommendedGap'],
                                     targetProbaRecommendedGap=kwargs['targetProbaRecommendedGap'],
                                     chromosomesRewrittenInTbs=arguments['mode:chromosomesRewrittenInTbs'],
                                     drawAllInformations=arguments['drawAllInformations'],
                                     scaleFactorRectangles=arguments['scaleFactorRectangles'],
                                     considerAllPairComps=arguments['considerAllPairComps'],
                                     switchOnDirectView=arguments['switchOnDirectView'],
                                     optimisation=kwargs['optimisation'],
                                     inSbsInPairComp=sbsInPairComp,
                                     outSyntenyBlocksFileName=arguments['out:SyntenyBlocks'],
                                     outImageFileName=arguments['out:ImageFileName'],
                                     verbose=arguments['verbose'])
