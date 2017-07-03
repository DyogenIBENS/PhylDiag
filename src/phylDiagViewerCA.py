#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr or hrc@ens.fr
# Licences GLP v3 and CeCILL v2

import utils.myLightGenomes as myLightGenomes
import utils.myTools as myTools
import utils.myGenomesDrawer as myGenomesDrawer

# import PhylDiag core algorithm
import utils.myDiags as myDiags

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
    [('withSbs', bool, True),
     ("in:syntenyBlocks", str, 'None'),
     ("out:syntenyBlocks", str, "./syntenyBlocksDrawer.txt"),
     ("mode:chromosomesRewrittenInTbs", bool, False),
     ('convertGenicToTbCoordinates', bool, False),
     ('drawAllInformations', bool, False),
     ("scaleFactorRectangles", float, 1.0),
     ("out:imageFileName", str, "./homologyMatrix.svg"),
     ("onlyROIcomp", bool, False),
     ('switchOnDirectView', bool, False),
     ('withIdsOfSbs', bool, True),
     ('verbose', bool, True)],
    __doc__)

kwargs = myDiags.defaultKwargsPhylDiag(arguments=arguments)
kwargs['verbose'] = True

genome1 = myLightGenomes.LightGenome(arguments['genome1'], withDict=True)
genome2 = myLightGenomes.LightGenome(arguments['genome2'], withDict=True)
families = myLightGenomes.Families(arguments['families'])

sbsInPairComp = None
if arguments['withSbs']:
    if arguments['in:syntenyBlocks'] != 'None':
        # load precomputed sbs if any
        sbsInPairComp = myDiags.parseSbsFile(arguments['in:syntenyBlocks'], genome1=genome1, genome2=genome2,
                                             withIds=arguments['withIdsOfSbs'])
        if arguments['withIdsOfSbs']:
            assert isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists)
        else:
            assert isinstance(sbsInPairComp, myTools.Dict2d)
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
                                     gapMax_Diag_Sbs=kwargs['gapMax_Diag_Sbs'],
                                     identifyMonoGenicCs=kwargs["identifyMonoGenicCs"],
                                     identifyMicroRearr=kwargs['identifyMicroRearr'],
                                     truncationMax=kwargs['truncationMax'],
                                     sameStrand=kwargs['sameStrand'],
                                     validateImpossToCalc_mThreshold=kwargs['validateImpossToCalc_mThreshold'],
                                     nbHpsRecommendedGap=kwargs['nbHpsRecommendedGap'],
                                     targetProbaRecommendedGap=kwargs['targetProbaRecommendedGap'],
                                     chromosomesRewrittenInTbs=arguments['mode:chromosomesRewrittenInTbs'],
                                     scaleFactorRectangles=arguments['scaleFactorRectangles'],
                                     onlyROIcomp=arguments['onlyROIcomp'],
                                     switchOnDirectView=arguments['switchOnDirectView'],
                                     optimisation=kwargs['optimisation'],
                                     inSbsInPairComp=sbsInPairComp,
                                     outSyntenyBlocksFileName=arguments['out:syntenyBlocks'],
                                     outImageFileName=arguments['out:imageFileName'],
                                     nbCaractersForGeneNamesInSymlbols=0,
                                     verbose=arguments['verbose'])

