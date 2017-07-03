#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr or hrc@ens.fr
# Licences GLP v3 and CeCILL v2

import os
import sys
import phylDiag
import argparse
import utils.myDiags as myDiags
import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes
import utils.myGenomesDrawer as myGenomesDrawer

p = argparse.ArgumentParser(description='Graphical visualisation of synteny blocks in homology matrices',
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                            parents=[phylDiag.p],
                            # FIXME
                            epilog='Warning: with --withoutSbs it may not return the desired homology matrix')
# positional arguments
p.add_argument('image', metavar='IMAGE', type=str, help='path to the returned image.svg')

# optional arguments
# if ROI1 and ROI2 are both None, it's a whole genome comparison, else its' a comparison of a pair of chromosomes
p.add_argument('--ROI1', metavar='chr1:beg1-end1', type=str, default=None, help='region of interest (ROI) on the first genome')
p.add_argument('--ROI2', metavar='chr2:beg2-end2', type=str, default=None, help='region of interest (ROI) on the second genome')

withSbs_parser = p.add_mutually_exclusive_group(required=False)
withSbs_parser.add_argument('--withSbs', dest='withSbs', action="store_true", help='draw sbs in the homology matrix')
withSbs_parser.add_argument('--withoutSbs', dest='withSbs', action="store_false")
p.set_defaults(withSbs=True)

p.add_argument('-s', '--inSbs', type=str, default=None, help='synteny blocks to draw (*.sbs)')
p.add_argument('-i', '--withSbIds', action="store_true", help='draw ids of synteny blocks')
p.add_argument('-b', '--geneIdxsToTbIdxs', action="store_true", help='convert gene idxs to tandem block idxs')
p.add_argument('-r', '--chrsInTbs', action="store_true", help='draw chromosomes in tandem blocks (after collapsing clusters of tandem duplicates)')
p.add_argument('-a', '--scaleRects', type=phylDiag.naturalFloat, default=1.0, help='scale factor of rectangle widths, if they are too small increase it')
# FIXME onlyROIcomp and ((ROI1 is None) and (ROI2 is None)) should not be authorized
# FIXME onlyROIcomp may not work properly
p.add_argument('--onlyROIcomp', action="store_true", help='consider only comparisons of both ROIs (change the filtering)')
p.add_argument('-l', '--liveView', action="store_true", help='turn on direct view with firefox as soon as the computation is finished')
p.add_argument('-o', '--outSbs', type=str, default='res/sbs.txt', help='information about drawn sbs')
args = p.parse_args()
myTools.printArguments(vars(args), sys.stderr)

if not args.imr and not args.imcs:
    print >> sys.stderr, 'no-imr and no-imcs: thus mmg is ignored'
if not args.truncation:
    print >> sys.stderr, 'no-truncation: thus truncationMax is ignored'

# FIXME, add a specific micro-gap for imr
gapMax_Diag_Sbs = phylDiag.setArgVal(args.imcs, args.mmg)
truncationMax = phylDiag.setArgVal(args.truncation, args.truncationMax)

genome1 = myLightGenomes.LightGenome(args.genome1, withDict=True)
genome2 = myLightGenomes.LightGenome(args.genome2, withDict=True)
families = myLightGenomes.Families(args.families)

sbsInPairComp = None
if args.onlyROIcomp:
    # sbs will be computed latter, when ROIs are compared
    pass
else:
    if args.withSbs:
        if args.inSbs is not None:
            # load precomputed sbs if any
            sbsInPairComp = myDiags.parseSbsFile(args.inSbs, genome1=genome1, genome2=genome2,
                                                 withIds=args.withSbIds)
            if args.withSbs:
                assert isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists)
            else:
                assert isinstance(sbsInPairComp, myTools.Dict2d)
        else:
            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, families,
                                                                filterType=args.filter,
                                                                tandemGapMax=args.tandemGapMax,
                                                                gapMax=args.gapMax,
                                                                distanceMetric=args.distanceMetric,
                                                                gapMax_Diag_Sbs=gapMax_Diag_Sbs,
                                                                identifyMonoGenicCs=args.imcs,
                                                                identifyMicroRearr=args.imr,
                                                                truncationMax=truncationMax,
                                                                verbose=args.verbose)
if (args.ROI1 is None) and (args.ROI2 is None):
    # FIXME with sbsInPairComp == None, this should return the matrix with only homologies
    try:
        os.remove(args.image)
    except:
        pass
    # whole genome comparison
    assert len(genome1.values()) > 0, "genome1 contains no chromosome, check the option \'removeUnofficialChromosomes\'"
    assert len(genome2.values()) > 0, "genome2 contains no chromosome, check the option \'removeUnofficialChromosomes\'"
    WHM = myGenomesDrawer.drawWholeGenomeHomologyMatrices(genome1, genome2, families,
                                                          # specific to phylDiag
                                                          filterType=args.filter,
                                                          minChromLength=args.minChrLen,
                                                          tandemGapMax=args.tandemGapMax,
                                                          # specific to the viewer
                                                          inSbsInPairComp=sbsInPairComp,
                                                          scaleFactorRectangles=args.scaleRects,
                                                          outputFileName=None,
                                                          maxWidth=100,
                                                          maxHeight=100,
                                                          # draw a black rectangle for each comp. where there is a sb
                                                          fillCompsWithSbs=False)
    if sbsInPairComp:
        nbSbs = len(sbsInPairComp.intoList())
    else:
        nbSbs = 0

    myGenomesDrawer.writeSVGFileForPairwiseCompOfGenomes(genome1.name,
                                                         genome2.name,
                                                         WHM,
                                                         # specific to phylDiag
                                                         filterType=args.filter,
                                                         tandemGapMax=args.tandemGapMax,
                                                         gapMax=args.gapMax,
                                                         distanceMetric=args.distanceMetric,
                                                         gapMax_Diag_Sbs=gapMax_Diag_Sbs,
                                                         identifyMicroRearr=args.imr,
                                                         truncationMax=truncationMax,
                                                         # specific to the viewer
                                                         chromosomesRewrittenInTbs=args.chrsInTbs,
                                                         nbSbs=nbSbs,
                                                         outImageFileName=args.image,
                                                         switchOnDirectView=args.liveView)
elif (args.ROI1 is not None) and (args.ROI2 is not None):
    # FIXME with sbsInPairComp == None, this should return the matrix with only homologies
    try:
        os.remove(args.image)
    except:
        pass
    # DEBUG
    # print >> sys.stderr, gapMax_Diag_Sbs
    # print >> sys.stderr, truncationMax
    # print >> sys.stderr, vars(args)
    # comparison of two ROIs
    myGenomesDrawer.homologyMatrixViewer(genome1, genome2, families, args.ROI1, args.ROI2,
                                         # specific to phylDiag
                                         filterType=args.filter,
                                         tandemGapMax=args.tandemGapMax,
                                         gapMax=args.gapMax,
                                         distanceMetric=args.distanceMetric,
                                         gapMax_Diag_Sbs=gapMax_Diag_Sbs,
                                         identifyMonoGenicCs=args.imcs,
                                         identifyMicroRearr=args.imr,
                                         truncationMax=truncationMax,
                                         verbose=args.verbose,
                                         # specific to the viewer
                                         convertGenicToTbCoordinates=args.geneIdxsToTbIdxs,
                                         minChromLength=args.minChrLen,
                                         chromosomesRewrittenInTbs=args.chrsInTbs,
                                         scaleFactorRectangles=args.scaleRects,
                                         onlyROIcomp=args.onlyROIcomp,
                                         switchOnDirectView=args.liveView,
                                         inSbsInPairComp=sbsInPairComp,
                                         outSyntenyBlocksFileName=args.outSbs,
                                         outImageFileName=args.image,
                                         nbCaractersForGeneNamesInSymlbols=0)
else:
    raise argparse.ArgumentError("One of the two ROI argument is not set in options")
