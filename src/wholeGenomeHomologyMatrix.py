#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# This script was originally written by Matthieu Muffato, it has been integrated
# to the PhylDiag package for practical uses

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
         ('outImageFileName', str, 'res/image.svg')],
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
                                                  scaleFactorRectangles=10,
                                                  outputFileName=None,
                                                  maxWidth=100,
                                                  maxHeight=100)

nbSbs = len(sbsInPairComp.intoList())

myGenomesDrawer.writeSVGFileForPairwiseCompOfGenomes(genome1.name,
                                                     genome2.name,
                                                     WHM,
                                                     chromosomesRewrittenInTbs=arguments['chromosomesRewrittenInTbs'],
                                                     filterType=kwargs['filterType'],
                                                     tandemGapMax=kwargs['tandemGapMax'],
                                                     gapMax=kwargs['gapMax'],
                                                     distanceMetric=kwargs['distanceMetric'],
                                                     gapMaxMicroInv=kwargs['gapMaxMicroInv'],
                                                     identifyBreakpointsWithinGaps=kwargs['gapMaxMicroInv'],
                                                     overlapMax=kwargs['overlapMax'],
                                                     nbSbs=nbSbs,
                                                     outImageFileName=arguments['outImageFileName'],
                                                     switchOnDirectView=False)