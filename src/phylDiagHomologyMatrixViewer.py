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

if __name__ == '__main__':

    arguments = myTools.checkArgs(
        [("genome1", file),
         ("genome2", file),
         ("families", file),
         ("chr1:deb1-fin1", str),
         ("chr2:deb2-fin2", str)],
        [("filterType", str, 'InBothGenomes'),
         ("minChromLength", int, 2),
         ("tandemGapMax", int, 0),
         ("distanceMetric", str, 'CD'),
         ("gapMax", str, 'None'),
         ("distinguishMonoGenicDiags", bool, True),
         ("pThreshold", str, 'None'),
         ('gapMaxMicroInv', str, '0'),
         ('identifyBreakpointsWithinGaps', bool, True),
         ("overlapMax", str, 'None'),
         ("consistentSwDType", bool, True),
         ("validateImpossToCalc_mThreshold", int, 3),
         ("in:SyntenyBlocks", str, 'None'),
         ("out:SyntenyBlocks", str, "./syntenyBlocksDrawer.txt"),
         ("mode:chromosomesRewrittenInTbs", bool, False),
         ('convertGenicToTbCoordinates', bool, False),
         ('nbHpsRecommendedGap', int, 2),
         ('targetProbaRecommendedGap', float, 0.01),
         ("scaleFactorRectangles", float, 2.0),
         ("out:ImageName", str, "./homologyMatrix.svg"),
         ("considerAllPairComps", bool, True),
         ('switchOnDirectView', bool, False),
         ('optimisation', str, 'cython'),
         ('verbose', bool, True)],
        __doc__)

for (argN, tpe) in [('gapMax', int), ('overlapMax', int), ('gapMaxMicroInv', int), ('pThreshold', float)]:
    if arguments[argN] == 'None':
        arguments[argN] = None
    else:
        try:
            arguments[argN] = tpe(arguments[argN])
        except:
            raise TypeError('%s is either an int or None' % argN)

gapMax = arguments['gapMax']
distinguishMonoGenicDiags = arguments['distinguishMonoGenicDiags']
gapMaxMicroInv = arguments['gapMaxMicroInv']
overlapMax = arguments['overlapMax']
tandemGapMax = arguments['tandemGapMax']
identifyBreakpointsWithinGaps = arguments['identifyBreakpointsWithinGaps']
consistentSwDType = arguments['consistentSwDType']
minChromLength = arguments['minChromLength']
pThreshold = arguments['pThreshold']
validateImpossToCalc_mThreshold = arguments['validateImpossToCalc_mThreshold']
chromosomesRewrittenInTbs = arguments['mode:chromosomesRewrittenInTbs']
convertGenicToTbCoordinates = arguments['convertGenicToTbCoordinates']
distanceMetric = arguments['distanceMetric']
nbHpsRecommendedGap = arguments['nbHpsRecommendedGap']
targetProbaRecommendedGap = arguments['targetProbaRecommendedGap']
considerAllPairComps = arguments['considerAllPairComps']
scaleFactorRectangles = arguments['scaleFactorRectangles']

# Load genomes
genome1 = myLightGenomes.LightGenome(arguments['genome1'], withDict=True)
genome2 = myLightGenomes.LightGenome(arguments['genome2'], withDict=True)
# load families
families = myLightGenomes.Families(arguments["families"])

inSbsInPairComp = None
if arguments['in:SyntenyBlocks'] == 'None':
    inSbsInPairComp = None
else:
    # load precomputed sbs if any
    inSbsInPairComp = arguments['in:SyntenyBlocks']
    inSbsInPairComp = myDiags.parseSbsFile(arguments['in:SyntenyBlocks'], genome1=genome1, genome2=genome2)

myGenomesDrawer.homologyMatrixViewer(genome1, genome2, families, arguments['chr1:deb1-fin1'], arguments['chr2:deb2-fin2'],
                                     convertGenicToTbCoordinates=arguments['convertGenicToTbCoordinates'],
                                     filterType=arguments['filterType'],
                                     minChromLength=arguments['minChromLength'],
                                     tandemGapMax=arguments['tandemGapMax'],
                                     distanceMetric=arguments['distanceMetric'],
                                     gapMax=arguments['gapMax'],
                                     distinguishMonoGenicDiags=arguments['distinguishMonoGenicDiags'],
                                     pThreshold=arguments['pThreshold'],
                                     gapMaxMicroInv=arguments['gapMaxMicroInv'],
                                     identifyBreakpointsWithinGaps=arguments['identifyBreakpointsWithinGaps'],
                                     overlapMax=arguments['overlapMax'],
                                     consistentSwDType=arguments['consistentSwDType'],
                                     validateImpossToCalc_mThreshold=arguments['validateImpossToCalc_mThreshold'],
                                     nbHpsRecommendedGap=arguments['nbHpsRecommendedGap'],
                                     targetProbaRecommendedGap=arguments['targetProbaRecommendedGap'],
                                     chromosomesRewrittenInTbs=arguments['mode:chromosomesRewrittenInTbs'],
                                     scaleFactorRectangles=arguments['scaleFactorRectangles'],
                                     considerAllPairComps=arguments['considerAllPairComps'],
                                     switchOnDirectView=arguments['switchOnDirectView'],
                                     optimisation=arguments['optimisation'],
                                     inSbsInPairComp=inSbsInPairComp,
                                     outSyntenyBlocksFileName=arguments['out:SyntenyBlocks'],
                                     outImageFileName=arguments['out:ImageName'],
                                     verbose=arguments['verbose'])