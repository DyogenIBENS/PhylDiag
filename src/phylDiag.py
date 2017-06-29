#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import argparse
import collections
import utils.myDiags as myDiags
import utils.myTools as myTools
import utils.myLightGenomes as myLightGenomes

def naturalInt(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError("Negative integer is not allowed")
    return x

def naturalFloat(x):
    x = float(x)
    if x < 0:
        raise argparse.ArgumentTypeError("Negative float are not allowed")
    return x

def setArgVal(argBool, argVal):
    if argBool:
        gapMax_Diag_Sbs = argVal
    else:
        gapMax_Diag_Sbs = None
        if args.mmg is not None:
            print >> sys.stderr, 'Warning imcs is False thus mmg is neglected for monog-genic conserved segments'
    return gapMax_Diag_Sbs

p = argparse.ArgumentParser(description='From the comparison of two extant genomes and corresponding gene families, PhylDiag detects conserved segments, i.e. segments of chromosomes unbroken during evolution.',
                            add_help=False,
                            )
# positional arguments
p.add_argument('genome1', metavar='G1', type=str, help='genome1')  # G1.genome.bz2
p.add_argument('genome2', metavar='G2', type=str, help='genome2')  # G2.genome.bz2'
p.add_argument('families', metavar='F', type=str, help='set of gene families')  # F.families.bz2

# optional arguments
## pre-processing
p.add_argument('-m', '--minChrLen', type=naturalInt, default=2, help='minimum number of genes in considered chromosomes')
p.add_argument('-f', '--filter', type=str, default=myDiags.FilterType.InBothGenomes, choices=list(myDiags.FilterType),
               help='filter type')
p.add_argument('-t', '--tandemGapMax', type=naturalInt, default=10,
               help='maximum gap between tandem duplicates in the same cluster')
p.add_argument('-d', '--distanceMetric', type=str, default=myDiags.DistanceMetric.CD, choices=myDiags.DistanceMetric,
               help='metric used for the calculation of 2D distances. CD: Chebyshev, MD: Manhattan, DPD: Diagonal Pseudo Distance, ED: Euclidian')
p.add_argument('-g', '--gapMax', type=naturalInt, default='5',
               help='maximum 2D gap between chained homologies')

## post-processing
imr_parser = p.add_mutually_exclusive_group(required=False)
imr_parser.add_argument('--imr', dest='imr', action='store_true',  help='identify micro-rearrangements')
imr_parser.add_argument('--no-imr', dest='imr', action='store_false')
p.set_defaults(imr=True)
imcs_parser = p.add_mutually_exclusive_group(required=False)
imcs_parser.add_argument('--imcs', dest='imcs', action='store_true',  help='identify mono-genic conserved segments')
imcs_parser.add_argument('--no-imcs', dest='imcs', action='store_false')
p.set_defaults(imcs=True)
p.add_argument('--mmg', type=naturalInt, default=1,
               help='maximum micro-gap, maximum gap allowed between: the homology of a detectable micro-segment and the nearest homology of a diagonal')
truncation_parser = p.add_mutually_exclusive_group(required=False)
truncation_parser.add_argument('--truncation', dest='truncation', action='store_true',  help='truncate overlapping diagonals')
truncation_parser.add_argument('--no-truncation', dest='truncation', action='store_false')
p.set_defaults(truncation=True)
p.add_argument('--truncationMax', type=naturalInt, default=10,
               help='maximum truncated length of the smallest overlapping diagonals, above the diag. is fully removed, without truncation')

## misc
### -v : set the verbosity level to 1
### -v X : set the verbosity level to X
### default is no verbosity
# p.add_argument("-v", "--verbosity", action="count", default=0,
#                help="increase output verbosity")
p.add_argument("-v", "--verbose", action="store_true", help="verbosity")

if __name__ == '__main__':
    pp = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[p])
    args = pp.parse_args()

    myTools.printArguments(vars(args), sys.stderr)

    if not args.imr and not args.imcs:
        print >> sys.stderr, 'no-imr and no-imcs: thus mmg is ignored'
    if not args.truncation:
        print >> sys.stderr, 'no-truncation: thus truncationMax is ignored'

    #identifyMicroRearr = True if args.imr else False
    #identifyMonoGenicCs = True if args.imcs else False
    gapMax_Diag_Sbs = setArgVal(args.imcs, args.mmg)
    truncationMax = setArgVal(args.truncation, args.truncationMax)

    genome1 = myLightGenomes.LightGenome(args.genome1)
    genome2 = myLightGenomes.LightGenome(args.genome2)
    families = myLightGenomes.Families(args.families)

    # extract synteny blocks
    sbsInPairComp = myDiags.extractSbsInPairCompGenomes(genome1, genome2, families,
                                                        minChromLength=args.minChrLen,
                                                        filterType=args.filter,
                                                        tandemGapMax=args.tandemGapMax,
                                                        gapMax=args.gapMax,
                                                        distanceMetric=args.distanceMetric,
                                                        gapMax_Diag_Sbs=gapMax_Diag_Sbs,
                                                        identifyMonoGenicCs=args.imcs,
                                                        identifyMicroRearr=args.imr,
                                                        truncationMax=truncationMax,
                                                        verbose=args.verbose)

    nbOfSbs = sum([len(sbsInPairComp[c1][c2]) for (c1, c2) in sbsInPairComp.keys2d()])
    print >> sys.stderr, "Number of synteny blocks = %s" % nbOfSbs

    if nbOfSbs > 0:
        sbLens = sorted([len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()])
        nbSbsOfLen = collections.Counter(sbLens)
        distribSbLens = [" %s:%s" % (nbSbsOfLen[sbLen], sbLen) for sbLen in sorted(nbSbsOfLen.keys())]
        distribSbLens = distribSbLens[:5] + ["..."] + distribSbLens[-3:]
        print >> sys.stderr, "Distribution of sb lengths (nbSbs:length) = %s" % ",".join(distribSbLens)

    # FIXME: create an empty file or do not create any file when no sbs
    myDiags.printSbsFile(sbsInPairComp, genome1, genome2, sortByDecrLengths=True)
