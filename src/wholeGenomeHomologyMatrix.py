#!/usr/bin/python
# -*- coding: utf-8 -*-

# PhylDiag version 2.0 (6/11/2015)
# python 2.7
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
        Plot the homology matrix between two genomes, with synteny blocks of PhylDiag
"""

from utils import myTools, myLightGenomes, myGenomesDrawer, myDiags

arguments = myTools.checkArgs(
        [("genome1", file), ("genome2", file), ("families", file)],
        myDiags.defaultArgsPhylDiag +\
        [('removeUnofficialChromosomes', bool, True),
         ('withSbs', bool, True),
         ('chromosomesRewrittenInTbs', bool, False),
         ('out:imageFileName', str, 'res/image.svg'),
         ('scaleFactorRectangles', float, 10.0)],
        __doc__)

kwargs = myDiags.defaultKwargsPhylDiag(arguments=arguments)
kwargs['verbose'] = True

genome1 = myLightGenomes.LightGenome(arguments['genome1'])
genome2 = myLightGenomes.LightGenome(arguments['genome2'])
families = myLightGenomes.Families(arguments['families'])

if arguments['removeUnofficialChromosomes']:
    genome1.removeUnofficialChromosomes()
    genome2.removeUnofficialChromosomes()

sbsInPairComp = None
if arguments['withSbs']:
    sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, families,
                                                        **kwargs)

WHM = myGenomesDrawer.drawWholeGenomeHomologyMatrices(genome1, genome2, families,
                                                  inSbsInPairComp=sbsInPairComp,
                                                  filterType=kwargs['filterType'],
                                                  minChromLength=kwargs['minChromLength'],
                                                  tandemGapMax=kwargs['tandemGapMax'],
                                                  scaleFactorRectangles=arguments['scaleFactorRectangles'],
                                                  outputFileName=None,
                                                  maxWidth=100,
                                                  maxHeight=100)

if sbsInPairComp:
    nbSbs = len(sbsInPairComp.intoList())
else:
    nbSbs = 0

myGenomesDrawer.writeSVGFileForPairwiseCompOfGenomes(genome1.name,
                                                     genome2.name,
                                                     WHM,
                                                     chromosomesRewrittenInTbs=arguments['chromosomesRewrittenInTbs'],
                                                     filterType=kwargs['filterType'],
                                                     tandemGapMax=kwargs['tandemGapMax'],
                                                     gapMax=kwargs['gapMax'],
                                                     distanceMetric=kwargs['distanceMetric'],
                                                     gapMaxMicroInv=kwargs['gapMaxMicroInv'],
                                                     identifyMicroRearrangements=kwargs['identifyMicroRearrangements'],
                                                     truncationMax=kwargs['truncationMax'],
                                                     nbSbs=nbSbs,
                                                     outImageFileName=arguments['out:imageFileName'],
                                                     switchOnDirectView=False)