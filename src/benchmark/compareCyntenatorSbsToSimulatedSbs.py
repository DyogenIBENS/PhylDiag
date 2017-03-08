#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import enum
import itertools
import pickle
import time
from utils import myTools
from utils import myIntervals
from utils import myLightGenomes
from utils import myDiags
from utils import myPhylTree
from utils import myMapping
from utils import myCyntenator
import collections

from libs import myBreakpointsAnalyser

# 3 line to change to adapt this script to a local installation
BIN_CYNTENATOR = '/home/jlucas/Libs/cyntenator/cyntenator'
DATA_CYNTENATOR = '/home/jlucas/Libs/PhylDiag/data/benchmark/cyntenator'
RES_CYNTENATOR = '/home/jlucas/Libs/PhylDiag/res/benchmark/cyntenator'

for directory in [DATA_CYNTENATOR, RES_CYNTENATOR]:
    if not os.path.exists(directory):
        os.makedirs(directory)

__doc__ = """This script compares the sbs of Cyntenator with the reference sbs recorded during simulation of magSimus
          (with the breakpoint analyser feature)
          """
FilterType = enum.Enum('InFamilies', 'InBothGenomes', 'None')

def analyseMistakeOfPhylDiag(geneNameInvolvedInAMistake, sbsInPairComp, simulatedSbsGenome,
                             genome1, genome2, mistake='Fp'):
    gn = geneNameInvolvedInAMistake
    assert mistake in {'Fp', 'Fn'}
    assert isinstance(simulatedSbsGenome, myLightGenomes.LightGenome)
    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    assert genome1.withDict is True
    assert genome2.withDict is True
    assert simulatedSbsGenome.withDict is True
    positionsOnGenomes = []
    if mistake == 'Fp':
        for ((c1, c2), sb) in sbsInPairComp.items2d():
            if geneNameInvolvedInAMistake == sb.la[0][0]:
                positionsOnGenomes.append((myLightGenomes.GeneP(c1, sb.l1[0][0]), myLightGenomes.GeneP(c2, sb.l2[0][0])))
            if geneNameInvolvedInAMistake == sb.la[-1][0]:
                positionsOnGenomes.append((myLightGenomes.GeneP(c1, sb.l1[-1][0]), myLightGenomes.GeneP(c2, sb.l2[-1][0])))
    elif mistake == 'Fn':
        posSbsGenome = simulatedSbsGenome.getPosition(gn, default=None)
        assert posSbsGenome is not None
        # le gene se situe a l'extremite d'un sb de reference (chromExtremity mode)
        assert posSbsGenome.idx in {0, len(simulatedSbsGenome[posSbsGenome.c]) - 1}
        positionsOnGenomes.append((genome1.getPosition(gn, default=None), genome2.getPosition(gn, default=None)))
    # print >> sys.stderr, "%s on G1:%s and G2:%s" % (mistake, str(possGenome1), str(possGenome2))
    return positionsOnGenomes

    # distribLenSbsWithFn = collections.defaultdict(int)
    # for (geneN, _) in sFn:
    #     geneP = sbsGenome.getPosition(geneN, default=None)
    #     if geneP:
    #         lenSbFn = len(sbsGenome[geneP.c])
    #         print >> sys.stderr, "Fp: gene %s in chrR %s of len=%s at idx %s" % (geneN, geneP.c, lenSbFn, geneP.idx)
    #         distribLenSbsWithFn[lenSbFn] += 1
    # sumLens = sum(distribLenSbsWithFn.values())
    # percentageLen1 = 100 * float(distribLenSbsWithFn[1]) / sumLens if sumLens != 0 else 'None'
    # print >> sys.stderr, "(%s%% of len1), distribLenSbsWithFn=%s" %\
    #                      (percentageLen1,
    #                       sorted([(l, nb) for (l, nb) in distribLenSbsWithFn.iteritems()], key=lambda x: x[0]))
    # return distribLenSbsWithFn

def editSbs(referenceSbs, sbs):
    if arguments['removeSbsOfLen1']:
        sbs.removeChrsStrictlySmallerThan(2)
    if arguments['removeSbsOfLen1']:
        assert all(len(sb) > 1 for sb in sbs.values())
        assert all(len(sb) > 1 for sb in referenceSbs.values())

    if arguments['reduceToSameGeneNames']:
        sgSbs = phylDiagSbsGenome.getGeneNames(checkNoDuplicates=False)
        sgReferenceSbs = referenceSbs.getGeneNames()
        # 1st reduce to the same set of genes
        setGenesInBoth = sgSbs & sgReferenceSbs
        setGenesToRemove = (sgSbs | sgReferenceSbs) - setGenesInBoth
        (sbs, _, (nbRemovedSbs1, nbRemovedGenes1)) =\
            myMapping.remapFilterGeneContent(phylDiagSbsGenome, setGenesToRemove)
        (nbRemovedSbs1, nbRemovedGenes1)
        (referenceSbs2, _, (nbRemovedSbs2, nbRemovedGenes2)) =\
             myMapping.remapFilterGeneContent(referenceSbs, setGenesToRemove)
        print >> sys.stderr, (nbRemovedSbs2, nbRemovedGenes2)
    else:
        referenceSbs2 = referenceSbs
    return (referenceSbs2, sbs)


arguments = myTools.checkArgs(
    [
        # 'data/speciesTree.phylTree'
        ('speciesTree', file),
        # 'Mus.musculus'
        ('species1', str),
        # 'Gallus.gallus'
        ('species2', str),
        # 'Amniota'
        ('ancGenes', str)
    ],
    [
        ('pSimGenomes', str, '../../res/simu1/genes.%s.list.bz2'),
        ('pAncGenes', str, '../../res/simu1/ancGenes.%s.list.bz2'),
        ('pSimulatedSbs', str, '../../res/simu1/sbs.genes.%s.%s.list.bz2'),
        # filterType should be in ['None', 'InFamilies', 'InBothGenomes']
        ('filterType', str, 'InBothGenomes'),
        ('modeOfComparison', str, 'chromExtremity'),

        #('gaps', str, '(-1,-2,-5,-20,-1000)'),
        ('gaps', str, '(-2,)'),
        ('mismatchs', str, '(-3,)'),
        # Les thresholds (-1 et 0) en dessous de 1 donnent le même résultat que le threshold 1
        # les valeurs intermédiaires 1.3, 1.7, ... semblent arrondies aux valeurs inférieures
        ('thresholds', str, '(1,2,3,4,5,10,50,100)'),
        #('thresholds', str, '(1, 1.3, 1.7, 2)'),
        ('coverage', int, 2),
        ('filter', int, 100),
        #('gapMaxs', str, '(1,5,20)'),
        #('gapMaxs', str, '(0,1)'),
        #('gapMaxs', str, '(0,)'),
        #('tandemGapMaxs', tuple, (0, 10, 100)),
        #('gapMaxMicroInv', int, 0),
        #('identifyMonoGenicInvs', bool, False),
        ('tandemGapMax', int, 5),
        ('distinguishMonoGenicDiags', bool, True),
        #('identifyMicroRearrangements', bool, False),

        #
        ('preComputePairwiseSbs', bool, True),
        #
        ('removeSbsOfLen1', bool, False),
        #
        ('reduceToSameGeneNames', bool, False),
        #
        ('oriented', bool, True),

        ('analyseMistakesOneByOne', bool, False),
        ('mistake', str, 'Fn'),
        ('outFigureName', str, './%s/Cyntenator.svg' % RES_CYNTENATOR)
    ],
    __doc__,
    # load Only Default Options
    loadOnlyDefaultOptions=False)
assert arguments['modeOfComparison'] in ['chromExtremity', 'adjacency']
# if arguments['gapMaxMicroInv'] == 'None':
#     arguments['gapMaxMicroInv'] = None
# else:
#     arguments['gapMaxMicroInv'] = int(arguments['gapMaxMicroInv'])
if arguments['gaps'] == 'None':
    arguments['gaps'] = None
else:
    arguments['gaps'] = tuple(eval(arguments['gaps']))
if arguments['mismatchs'] == 'None':
    arguments['mismatchs'] = None
else:
    arguments['mismatchs'] = tuple(eval(arguments['mismatchs']))
if arguments['thresholds'] == 'None':
    arguments['thresholds'] = None
else:
    arguments['thresholds'] = tuple(eval(arguments['thresholds']))

# import os
# os.chdir('/home/jlucas/Libs/MagSimus')
# arguments['speciesTree'] = 'data/speciesTree.phylTree'
# arguments['species1'] = 'Mus.musculus'
# arguments['species2'] = 'Gallus.gallus'
# arguments['ancGenes'] = 'Amniota'
# arguments['filterType'] = 'InFamilies'
# arguments['preComputePairwiseSbs'] = False
# arguments['tandemGapMax'] = 0
# arguments['removeSbsOfLen1'] = True
# arguments['reduceToSameGeneNames'] = False

speciesTree = myPhylTree.PhylogeneticTree(arguments['speciesTree'])

if arguments['preComputePairwiseSbs']:
    myBreakpointsAnalyser.computePairwiseSyntenyBlocksOfMagSimus(speciesTree, arguments['pSimGenomes'],
                                                                 arguments['pAncGenes'], arguments['pSimulatedSbs'])

genome1 = myLightGenomes.LightGenome(arguments['pSimGenomes'] % arguments['species1'])
genome2 = myLightGenomes.LightGenome(arguments['pSimGenomes'] % arguments['species2'])
ancGenes = myLightGenomes.Families(arguments['pAncGenes'] % arguments['ancGenes'])
filterType = list(myDiags.FilterType._keys)
filterType = myDiags.FilterType[filterType.index(arguments["filterType"])]
# to have the species combination in alphabetical order
speciesNames = sorted([arguments['species1'], arguments['species2']])
referenceSbs = myLightGenomes.LightGenome(arguments['pSimulatedSbs'] % (speciesNames[0], speciesNames[1]))

if arguments['removeSbsOfLen1']:
    referenceSbs.removeChrsStrictlySmallerThan(2)
nbSbsR = len(referenceSbs.keys())
dictStatsCompPhylDiag = {}
dictStatsCompAdhore = {}
#pThresholdPhylDiag = arguments['pThresholdPhylDiag'] if arguments['pThresholdPhylDiag'] != 'None' else None


variableArg = ('ibwg-imgi-gmmi-tm', [(False, False, 0, None),
                                     (True, False, 0, None),
                                     (False, True, 0, None),
                                     (True, True, 0, None),
                                     (True, True, 1, None),
                                    (True, True, 1, 4)])

# variableArg = ('ibwg-imgi-gmmi-tm',
# [(True, True, 0, None),
# (True, True, 1, None)])

#variableArg = ('ibwg-imgi-tm', [(False, False, None)])
# tandemGapMaxs = arguments['tandemGapMaxs']
# identifyMonoGenicInversionS = [False, True]
# pThresholdS = [1.0, 0.1, 0.01, 0.001, None]
# gapMaxMicroInvS = [0, 1, 2, 3, 4, 5]
# distanceMetricS = ['CD', 'MD', 'DPD', 'ED']

# For the analysis of Fp and Fn
# ((g1_tb, mtb2g1, (nCL1, nGL1)), (g2_tb, mtb2g2, (nCL2, nGL2))) =\
#     myDiags.editGenomes(genome1, genome2, ancGenes,
#                         filterType=filterType, labelWith='FamName',
#                         tandemGapMax=arguments['tandemGapMax'], keepOriginal=True)
# g1_tb.computeDictG2P()
# g1_tb.computeDictG2Ps()
# g2_tb.computeDictG2P()
# g2_tb.computeDictG2Ps()
#######

mistake = arguments['mistake']

# be sure that the names differ between the two species
for chr in genome1.keys():
    genome1[chr] = [myLightGenomes.OGene('1'+g.n, g.s) for g in genome1[chr]]
for chr in genome2.keys():
    genome2[chr] = [myLightGenomes.OGene('2'+g.n, g.s) for g in genome2[chr]]
newFamilies = myLightGenomes.Families()
for family in ancGenes:
    newDns = set()
    for dn in family.dns:
        newDns.add('1'+dn)
        newDns.add('2'+dn)
    newFamilies.addFamily(myLightGenomes.Family(family.fn, newDns))

progressBar = myTools.ProgressBar(len(arguments['gaps']) * len(arguments['mismatchs']) * len(arguments['thresholds']) )
currLen = 1
dictStatsCompCyntenator = collections.defaultdict(tuple)
resultPickled = []
for gap in arguments['gaps']:
    print >> sys.stderr, "gap = %s" % gap
    for mismatch in arguments['mismatchs']:
        print >> sys.stderr, "mismatch = %s" % mismatch
        for threshold in arguments['thresholds']:
            print >> sys.stderr, "threshold = %s" % threshold
            currLen += 1
            progressBar.printProgressIn(sys.stderr, currLen, prefixLabel='comp. cyn. ')
            cyntenatorSbsGenome = myCyntenator.launchCyntenator(genome1, genome2, newFamilies,
                             # options of cyntenator
                             threshold=threshold, gap=gap, mismatch=mismatch,
                             coverage=arguments['coverage'], filter=arguments['filter'],
                             # options more general
                             tandemGapMax=arguments['tandemGapMax'],
                             filterType=myDiags.FilterType.InBothGenomes,
                             outCyntenatorGenomes="%s/genome.%s.list.cyntenator" % (DATA_CYNTENATOR, '%s'),
                             pathCyntenatorBin=BIN_CYNTENATOR,
                             outAlignmentCyntenator='%s/alignment.cyntenator' % RES_CYNTENATOR)

            (referenceSbs2, cyntenatorSbsGenome) = editSbs(referenceSbs, cyntenatorSbsGenome)

            nbCyntenatorSbs = len(cyntenatorSbsGenome.keys())
            try:
                (eff, (sTp, sFp, sFn)) = myIntervals.compareGenomes(cyntenatorSbsGenome, referenceSbs2,
                                                                    mode=arguments['modeOfComparison'], oriented=arguments['oriented'], verbose=True)
                'chromExtremity'
                sTpNs = set(gn for (gn, _) in sTp) if arguments['oriented'] else sTp
                dictStatsCompCyntenator[(gap, mismatch, threshold)] = (eff, (nbCyntenatorSbs, nbSbsR))
                # serialize results
                resultPickled.append( ((gap, mismatch, threshold), (eff.sn, eff.sp), (nbCyntenatorSbs, nbSbsR)) )
                tic = time.time()
                pickle.dump(resultPickled, open('%s/cyntenator_param%s.txt' % (RES_CYNTENATOR, currLen), 'w'))
                tac = time.time()
                print >> sys.stderr, 'resultPickled has been serialized in cyntenator_param%s.txt (it took %s s)' % (currLen, tac - tic)
                print >> sys.stderr, dictStatsCompCyntenator
            except:
                print >> sys.stderr, 'one computation was impossible !!!'
                pass

# # deserialize results
#currLen = 1
#dictStatsCompCyntenator = pickle.load(open('%s/cyntenator_param%s.txt' % (RES_CYNTENATOR, currLen), 'r'))

print dictStatsCompCyntenator
dictStatsCompCyntenator = dict([(k,(effs,foo)) for (k,effs,foo) in dictStatsCompCyntenator])
dictStatsCompCyntenator = dict([(k,v) for (k,v) in dictStatsCompCyntenator.iteritems() if k[2] in arguments['thresholds']])

####################################################################
# Plot informations
####################################################################

import matplotlib.pyplot as plt
def plot(dictStatsComp, arguments, variableArg):

    algo = 'Cyntenator'
    label = algo
    g = arguments['gaps'][0]
    m= arguments['mismatchs'][0]

    fig, axes = plt.subplots(nrows=2, ncols=2)
    axis_sensitivity = axes[0][0]
    axis_specificity = axes[1][0]
    axis_nbConservedSegments = axes[1][1]
    axis_legend = axes[0][1]
    thresholds = arguments['thresholds']
    # fig, axes = plt.subplots(nrows=5)
    # Xs
    # lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linecycler = itertools.cycle(lines)
    colorcycler = itertools.cycle(colors)
    assert len(arguments['gaps']) == 1 and len(arguments['mismatchs']) == 1
    supTitle = "Efficiency analysis %s - %s\n" % (arguments['species1'], arguments['species2']) +\
               'filterType=%s\n' % arguments['filterType'] + \
               'tandemGapMax=%s\n' % arguments['tandemGapMax'] + \
               'gaps=%s\n' % g +\
               'mismatchs=%s\n' % m
               #'distanceMetric=%s\n' % arguments['distanceMetric'] +\
               #'identifyMicroRearrangements=%s\n' % arguments['identifyMicroRearrangements'] +\
               #'removeSbsOfLen1=%s\n' % arguments['removeSbsOfLen1'] +\
               #'pThreshold=%s\n' % arguments['pThreshold'] +\
               #'gapMaxMicroInv=%s\n' % arguments['gapMaxMicroInv'] +\
               #'identifyMonoGenicInvs=%s\n' % arguments['identifyMonoGenicInvs'] +\

    fig.suptitle(supTitle, fontsize=12, )
    sensitivityTitle = 'sensitivity'
    specificityTitle = 'specificity'
    subTitles = (sensitivityTitle, specificityTitle)
    minY = {sensitivityTitle: sys.maxint, specificityTitle: sys.maxint}
    print >> sys.stderr, dictStatsComp

    color = next(colorcycler)
    line = next(linecycler)
    statsPerGapMax = dictStatsComp
    sensitivitys = [sn for ((sn, sp), (nbSbs, nbSbsR )) in [statsPerGapMax[(g,m,t)] for t in thresholds]]
    specificitys = [sp for ((sn, sp), (nbSbs, nbSbsR )) in [statsPerGapMax[(g,m,t)] for t in thresholds]]
    Yss = (sensitivitys, specificitys)
    for subTitle, ax, Ys in zip(subTitles, (axis_sensitivity, axis_specificity), Yss):
        minY[subTitle] = min(Ys) if min(Ys) < minY[subTitle] else minY[subTitle]
        print len(Ys)
        print len(thresholds)
        ax.plot(xrange(len(thresholds)), Ys, color=color, linestyle=line, label=label)
        ax.set_ylabel(subTitle)
        #ax.set_ylim((minY[subTitle], 1.0))
        ax.set_ylim((0.0, 1.0))
        ax.set_xticks(xrange(len(thresholds)))
        ax.set_xticklabels([str(threshold) for threshold in thresholds])
        ax.set_xlabel('thresholds')
    nbSbss = [nbSbs for (eff, (nbSbs, nbSbsR )) in [statsPerGapMax[(g,m,t)] for t in thresholds]]
    axis_nbConservedSegments.plot(xrange(len(thresholds)), nbSbss, color=color, linestyle=line, label=label)
    axis_nbConservedSegments.set_ylabel('#conserved segments')
    axis_nbConservedSegments.set_xlabel('thresholds')
    axis_nbConservedSegments.set_xticks(xrange(len(thresholds)))
    axis_nbConservedSegments.set_xticklabels([str(threshold) for threshold in thresholds])

    nbSbsEs = [nbSbsR for (eff, (nbSbs, nbSbsR)) in [statsPerGapMax[(g,m,t)] for t in thresholds]]
    #line2, = axes[-1].plot(gapMaxs, nbSbsEs, color='red', linestyle='-', label='Sim', linewidth=2)
    line2, = axis_nbConservedSegments.plot(xrange(len(thresholds)), nbSbsEs, color='red', linestyle='-', linewidth=1.5)
    legend2 = axis_nbConservedSegments.legend([line2], ('Simulation',), loc="upper right")
    #leg = axis_legend.legend(loc='center', ncol=2, shadow=True, title=variableArg[0], fancybox=True)
    leg = axis_nbConservedSegments.legend(loc='center', shadow=True, fancybox=True)
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.gca().add_artist(legend2)
    fig.tight_layout()
    plt.show(block=True)
    plt.savefig(arguments['outFigureName'], format='svg')

plot(dictStatsCompCyntenator, arguments, variableArg)