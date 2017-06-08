#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys

import itertools

import numpy

import enum
import time
import pickle
import getpass # https://stackoverflow.com/questions/842059/is-there-a-portable-way-to-get-the-current-username-in-python
import collections
from utils import myTools
from utils import myDiags
from utils import myPhylTree
from utils import myIntervals
from utils import myCyntenator
from utils import myLightGenomes

# TODO When MagSimus is published uncomment the next line, remove benchTools.computePairwiseSyntenyBlocksOfMagSimus,
# remove its sub functions and update links properly.
# from libs import myBreakpointsAnalyser
import benchmarkTools as benchTools

LOAD_PRECOMPUTED_BENCHMARK = True
# 5 lines to change to adapt this script to a local installation
BIN_CYNTENATOR = myCyntenator.PATH_CYNTENATOR_BIN

#DATA_CYNTENATOR = '/home/jlucas/Libs/PhylDiag/data/benchmark/cyntenator'
# You might want to change this folder depending on where you installed PhylDiag
# /home/<user>/Libs/PhylDiag/res/benchmark/cyntenator
RES_CYNTENATOR = '/home/' + getpass.getuser() + '/Libs/PhylDiag/res/benchmark/cyntenator'
DATA_CYNTENATOR = RES_CYNTENATOR
serialisationFile = RES_CYNTENATOR + '/serialisationCyntenatorBenchmarkValues'


for directory in [DATA_CYNTENATOR, RES_CYNTENATOR]:
    if not os.path.exists(directory):
        os.makedirs(directory)

__doc__ = """This script compares the RSBs of Cyntenator with the reference conserved segments (cs) recorded during
the simulation of MagSimus (with the breakpoint analyser feature)
"""
FilterType = enum.Enum('InFamilies', 'InBothGenomes', 'None')

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
        # data (in PhylDiag/data/benchmark)
        ('pSimGenomes', str, 'genes.%s.list.bz2'),
        ('pAncGenes', str, 'ancGenes.%s.list.bz2'),
        ('pSimulatedSbs', str, 'cs.genes.%s.%s.list.bz2'),

        # pre-process
        # filterType should be in ['None', 'InFamilies', 'InBothGenomes']
        ('filterType', str, 'InBothGenomes'),
        ('tandemGapMax', int, 5),

        # mode of comparison : 'adjacency' or 'chromExtremity' or 'familyContent'
        ('oriented', bool, True),

        # parameters (may contain several values, varying)
        #('gaps', str, '(-1,-2,-5,-20,-1000)'),
        #('gaps', str, '(-2,)'),
        ('gaps', str, '(-2, -1000000000)'),
        ('mismatchs', str, '(-3, -100000000)'),
        ('coverage', int, 2),
        ('filter', int, 100),
        #('thresholds', str, '(1,)'),
        ('thresholds', str, '(1,2,3,5,10,50)'),
        #('thresholds', str, '(1,10)'),
        # Les thresholds (..., -1000, ..., -1 et 0) donnent le même résultat que le threshold 1 !
        # les valeurs intermédiaires 1.3, 1.7, ... semblent arrondies aux valeurs inférieures
        #('thresholds', str, '(1, 1.3, 1.7, 2)'),
        ('includeMicroRearrangedGenesInSbs', bool, True),

        ('preComputePairwiseSbs', bool, True),

        # edit sbs found
        ('removeSbsOfLen1', bool, False),
        ('reduceToSameGeneNames', bool, False),

        #('analyseMistakesOneByOne', bool, False),
        #('mistake', str, 'Fn'),

        ('outFigureName', str, '%s/Cyntenator.svg' % RES_CYNTENATOR)
    ],
    __doc__,
    # load Only Default Options
    loadOnlyDefaultOptions=False)

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

if not LOAD_PRECOMPUTED_BENCHMARK :
    filterType = list(myDiags.FilterType._keys)
    filterType = myDiags.FilterType[filterType.index(arguments["filterType"])]

    # Load data
    speciesTree = myPhylTree.PhylogeneticTree(arguments['speciesTree'])
    if arguments['preComputePairwiseSbs']:
        benchTools.computePairwiseSyntenyBlocksOfMagSimus(speciesTree, arguments['pSimGenomes'],
                                                                     arguments['pAncGenes'], arguments['pSimulatedSbs'])
    genome1 = myLightGenomes.LightGenome(arguments['pSimGenomes'] % arguments['species1'])
    genome2 = myLightGenomes.LightGenome(arguments['pSimGenomes'] % arguments['species2'])
    ancGenes = myLightGenomes.Families(arguments['pAncGenes'] % arguments['ancGenes'])
    # to have the species combination in alphabetical order
    speciesNames = sorted([arguments['species1'], arguments['species2']])
    sbsTrue = myLightGenomes.LightGenome(arguments['pSimulatedSbs'] % (speciesNames[0], speciesNames[1]))

    # edit ref sbs
    if arguments['removeSbsOfLen1']:
        sbsTrue.removeChrsStrictlySmallerThan(2)
    nbSbsTrue = len(sbsTrue.keys())
    nbAncGenesInSbsTrue = len(sbsTrue.getGeneNames(asA=set, checkNoDuplicates=True))

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

    dictStatsCompPhylDiag = {}
    dictStatsCompAdhore = {}

    (newGenome1, newGenome2, newFamilies) = benchTools.ensureNamesDifferBetweenGenome1AndGenome2(genome1, genome2, ancGenes)

    progressBar = myTools.ProgressBar(len(arguments['gaps']) * len(arguments['mismatchs']) * len(arguments['thresholds']) )
    currLen = 1
    dictStatsCompCyntenator = {}
    resultPickled = []
    for gap in arguments['gaps']:
        print >> sys.stderr, "gap = %s" % gap
        for mismatch in arguments['mismatchs']:
            print >> sys.stderr, "mismatch = %s" % mismatch
            for threshold in arguments['thresholds']:
                suffix = 'g%s_m%s_t%s' %  (gap, mismatch, threshold)
                dictStatsCompCyntenator[(gap, mismatch, threshold)] = {}
                print >> sys.stderr, "threshold = %s" % threshold
                currLen += 1
                progressBar.printProgressIn(sys.stderr, currLen, prefixLabel='comp. cyn. ')
                cyntenatorSbsGenome = myCyntenator.launchCyntenator(newGenome1, newGenome2, newFamilies,
                                                                    # options of cyntenator
                                                                    threshold=threshold, gap=gap, mismatch=mismatch,
                                                                    coverage=arguments['coverage'],
                                                                    filter=arguments['filter'],
                                                                    # options more general
                                                                    tandemGapMax=arguments['tandemGapMax'],
                                                                    filterType=myDiags.FilterType.InBothGenomes,
                                                                    includeMicroRearrangedGenesInSbs=arguments['includeMicroRearrangedGenesInSbs'],
                                                                    outCyntenatorGenomes="%s/genome.%s.list.cyntenator" % (
                                                                    DATA_CYNTENATOR, '%s'),
                                                                    pathCyntenatorBin=BIN_CYNTENATOR,
                                                                    outAlignmentCyntenator='%s/alignment.cyntenator_%s' % (RES_CYNTENATOR, suffix))

                (sbsTrue1, cyntenatorSbsGenome) = benchTools.editSbs(sbsTrue, cyntenatorSbsGenome,
                                                                     removeSbsOfLen1=arguments['removeSbsOfLen1'],
                                                                     reduceToSameGeneNames=arguments['reduceToSameGeneNames'])

                nbCyntenatorSbs = len(cyntenatorSbsGenome.keys())
                nbFamiliesInSbsCyntenator = len(set(cyntenatorSbsGenome.getGeneNames(asA=set, checkNoDuplicates=False)))
                print >> sys.stderr, "analyse ancestral genes content of cyntenator sbs"
                benchTools.analyseFamiliesInSyntenyBlocks(cyntenatorSbsGenome, sbsTrue1)

                comp_adjacency = myIntervals.compareGenomes(cyntenatorSbsGenome, sbsTrue1,
                                                            mode='adjacency',
                                                            oriented=arguments['oriented'], verbose=True, returnSets=True)
                comp_chromExtremity = myIntervals.compareGenomes(cyntenatorSbsGenome, sbsTrue1,
                                                                 mode='chromExtremity',
                                                                 oriented=arguments['oriented'], verbose=True, returnSets=True)
                comp_familyContent = myIntervals.compareGenomes(cyntenatorSbsGenome, sbsTrue1,
                                                              mode='geneName',
                                                              oriented=arguments['oriented'], verbose=True, returnSets=True)
                dictStatsCompCyntenator[(gap, mismatch, threshold)]['adjacency'] = comp_adjacency[0]
                dictStatsCompCyntenator[(gap, mismatch, threshold)]['chromExtremity'] = comp_chromExtremity[0]
                dictStatsCompCyntenator[(gap, mismatch, threshold)]['familyContent'] = comp_familyContent[0]
                # do not serialize this (too huge)
                # dictStatsCompCyntenator[(gap, mismatch, threshold)]['adjacency_sets'] = comp_adjacency[1]
                # dictStatsCompCyntenator[(gap, mismatch, threshold)]['chromExtremity_sets'] = comp_chromExtremity[1]
                # dictStatsCompCyntenator[(gap, mismatch, threshold)]['familyContent_sets'] = comp_familyContent[1]
                dictStatsCompCyntenator[(gap, mismatch, threshold)]['nbSbs'] = (nbCyntenatorSbs, nbSbsTrue)
                dictStatsCompCyntenator[(gap, mismatch, threshold)]['nbFamilies'] = (nbFamiliesInSbsCyntenator, nbAncGenesInSbsTrue)

                # serialize results (inside the for loops because cyntenator is slow)
                cyntenatorSbsGenome.printIn(open(RES_CYNTENATOR + '/genome.sbs_%s.list.bz2' % suffix, 'w'))
                if os.path.isfile(serialisationFile):
                    os.remove(serialisationFile)
                pickle.dump(dictStatsCompCyntenator, open(serialisationFile, 'w'))
else:
    # deserialize results
    dictStatsCompCyntenator = pickle.load(open(serialisationFile, 'r'))

####################################################################
# Plot Results
####################################################################
import matplotlib.pyplot as plt
def plotRecallAndPrecision(dictStatsComp, arguments):
    fig, axes = plt.subplots(nrows=2, ncols=3)
    colTitles = ['cs extremities', 'cs genes adjacencies', 'cs anc. families']
    rowTitles = ['sensitivity', 'specificity']

    # varying parameter
    thresholds = arguments['thresholds']
    #thresholds = (1,2,3,4,5,10,50)

    for ax, colTitle in zip(axes[0], colTitles):
        ax.set_title(colTitle)
    for ax, rowTitle in zip(axes[:, 0], rowTitles):
        ax.set_ylabel(rowTitle, rotation=90, size='large')
    for ax in axes[0]:
        ax.set_xticks(xrange(len(thresholds)))
        ax.set_xticklabels([str(threshold) for threshold in thresholds])
    for ax in axes[1]:
        ax.set_xticks(xrange(len(thresholds)))
        ax.set_xticklabels([str(threshold) for threshold in thresholds])
        ax.set_xlabel('threshold')
    lineWidth = 2.5
    # for (i,j), ax in numpy.ndenumerate(axes):
    #     ax.set_ylim((0, 1))
    axes[0][0].set_ylim((0, 1))
    axes[1][0].set_ylim((0, 1))

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    colorcycler = itertools.cycle(colors)
    # lines = ["-", "--", "-.", ":"]
    # linecycler = itertools.cycle(lines)
    lineStyle = '-'

    # choose values of parameters to plot
    m = -3
    linesCyntenator= [None] * (len(arguments['gaps'])* len(arguments['mismatchs']))
    cpt = 0
    for ig, g in enumerate(arguments['gaps']):
        for im, m in enumerate(arguments['mismatchs']):
            assert g in arguments['gaps']
            assert m in arguments['mismatchs']

            color = next(colorcycler)
            label = 'Cyntenator (%s, %s)' % (g, m)

            modeOfComparisons = ['chromExtremity', 'adjacency', 'familyContent']
            for i, modeOfComparison in enumerate(modeOfComparisons):
                axis_sensitivity = axes[0][i]
                axis_specificity = axes[1][i]
                sensitivitys = [eff.sn for eff in [dictStatsComp[(g,m,t)][modeOfComparison] for t in thresholds]]
                specificitys = [eff.sp for eff in [dictStatsComp[(g,m,t)][modeOfComparison] for t in thresholds]]
                Yss = (sensitivitys, specificitys)
                for ax, Ys in zip((axis_sensitivity, axis_specificity), Yss):
                    linesCyntenator[cpt], = ax.plot(xrange(len(thresholds)), Ys, color=color, linestyle=lineStyle, linewidth=lineWidth)
            cpt+=1
    # nbSbss = [nbSbs for (_, (nbSbs, _), _) in [dictStatsComp[(g,m,t)] for t in thresholds]]
    # nbSbsTrues = [nbSbsTrue for (eff, (_, nbSbsTrue), _) in [dictStatsComp[(g, m, t)] for t in thresholds]]
    # nbFamiliesInSbss = [nbFamiliesInSbs for (_, _, (nbFamiliesInSbs, _)) in [dictStatsComp[(g, m, t)] for t in thresholds]]
    # nbAncGenesInSbsTrues = [nbAncGenesInSbsTrue for (_, _, (_, nbAncGenesInSbsTrue)) in [dictStatsComp[(g, m, t)] for t in thresholds]]
    #                                        linewidth=3.0)

    # change the format of y values
    for (i, j), ax in numpy.ndenumerate(axes):
        vals = ax.get_yticks()
        if (i == 1 and j == 2):
            ax.set_yticklabels(['{:3.1f}%'.format(x * 100) for x in vals])
        else:
            ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    lines = linesCyntenator
    labels = ['Cyntenator (%s, %s)' % (g, m) for g in arguments['gaps'] for m in arguments['mismatchs']]
    titleLegend = '(gaps, mismatchs)'
    assert len(lines) == len(labels), '%s = %s' % (len(lines), len(labels))
    fig.legend(lines, labels, ncol=2, title=titleLegend, loc='upper center', fontsize=15)

    fig.tight_layout()
    plt.show(block=True)
    plt.savefig(arguments['outFigureName'], format='svg')

def plotRecallPrecisionAndF1(dictStatsComp, arguments, Mylimits):
    fig, axes = plt.subplots(nrows=3, ncols=3)
    colTitles = ['extremities', 'adjacencies', 'gene names']
    rowTitles = ['recall', 'precision', 'F1-score']

    # varying parameter
    thresholds = arguments['thresholds']
    # thresholds = (1,2,3,4,5,10,50)

    for ax, colTitle in zip(axes[0], colTitles):
        ax.set_title(colTitle)
    for ax, rowTitle in zip(axes[:, 0], rowTitles):
        ax.set_ylabel(rowTitle, rotation=90, size='large')
    for ax in axes.flatten():
        ax.set_xticks(xrange(len(thresholds)))
        ax.set_xticklabels([str(threshold) for threshold in thresholds])
    for ax in axes[2]:
        ax.set_xlabel('threshold')
    lineWidth = 2.5
    # for (i,j), ax in numpy.ndenumerate(axes):
    #     ax.set_ylim((0, 1))
    # Axes limits
    # axes[0][0].set_ylim((0, 1))
    # axes[1][0].set_ylim((0, 1))
    # axes[2][0].set_ylim((0, 1))
    # axes[2][2].set_ylim((0.92, 1))
    # axes[0][2].set_ylim((0.95, 1))
    for (i, j), yl in numpy.ndenumerate(Mylimits):
        axes[i][j].set_ylim(yl)

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    colorcycler = itertools.cycle(colors)
    # lines = ["-", "--", "-.", ":"]
    # linecycler = itertools.cycle(lines)
    lineStyle = '-'

    # choose values of parameters to plot
    m = -3
    linesCyntenator = [None] * (len(arguments['gaps']) * len(arguments['mismatchs']))
    cpt = 0
    for ig, g in enumerate(arguments['gaps']):
        for im, m in enumerate(arguments['mismatchs']):
            assert g in arguments['gaps']
            assert m in arguments['mismatchs']

            color = next(colorcycler)
            label = 'Cyntenator (%s, %s)' % (g, m)

            modeOfComparisons = ['chromExtremity', 'adjacency', 'familyContent']
            for i, modeOfComparison in enumerate(modeOfComparisons):
                axis_recall = axes[0][i]
                axis_precision = axes[1][i]
                axis_f1 = axes[2][i]
                recall = [eff.r for eff in [dictStatsComp[(g, m, t)][modeOfComparison] for t in thresholds]]
                precision = [eff.p for eff in [dictStatsComp[(g, m, t)][modeOfComparison] for t in thresholds]]
                f1s = [eff.f1 for eff in [dictStatsComp[(g, m, t)][modeOfComparison] for t in thresholds]]
                Yss = (recall, precision, f1s)
                for ax, Ys in zip((axis_recall, axis_precision, axis_f1), Yss):
                    linesCyntenator[cpt], = ax.plot(xrange(len(thresholds)), Ys, color=color, linestyle=lineStyle,
                                                    linewidth=lineWidth)
            cpt += 1
    # nbSbss = [nbSbs for (_, (nbSbs, _), _) in [dictStatsComp[(g,m,t)] for t in thresholds]]
    # nbSbsTrues = [nbSbsTrue for (eff, (_, nbSbsTrue), _) in [dictStatsComp[(g, m, t)] for t in thresholds]]
    # nbFamiliesInSbss = [nbFamiliesInSbs for (_, _, (nbFamiliesInSbs, _)) in [dictStatsComp[(g, m, t)] for t in thresholds]]
    # nbAncGenesInSbsTrues = [nbAncGenesInSbsTrue for (_, _, (_, nbAncGenesInSbsTrue)) in [dictStatsComp[(g, m, t)] for t in thresholds]]
    #                                        linewidth=3.0)

    # change the format of x and y values
    #setOfAxesCoordsWithIncreasedPrecision = {(1, 2)}
    setOfAxesCoordsWithIncreasedPrecision = {}
    for (i, j), ax in numpy.ndenumerate(axes):
        ax.tick_params(labelsize=10)
        vals = ax.get_yticks()
        if not (i, j) in setOfAxesCoordsWithIncreasedPrecision:
            ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])
        else:
            ax.set_yticklabels(['{:3.1f}%'.format(x * 100) for x in vals])

    lines = linesCyntenator
    labels = ['Cyntenator (%s, %s)' % (g, m) for g in arguments['gaps'] for m in arguments['mismatchs']]
    titleLegend = '(gaps, mismatchs)'
    assert len(lines) == len(labels), '%s = %s' % (len(lines), len(labels))
    fig.legend(lines, labels, ncol=2, title=titleLegend, loc='upper center', fontsize=12)

    # print ylim
    for (i, j), ax in numpy.ndenumerate(axes):
        (yl, yh) = ax.get_ylim()
        print >> sys.stderr, "ax[%s,%s].ylim = (%.2f, %.2f)" % (i,j,yl,yh)

    fig.tight_layout()
    plt.show(block=True)
    plt.savefig(arguments['outFigureName'], format='svg')

arguments['mismatchs'] = [-3]
arguments['thresholds'] = (1,2,3,5,10)
plotRecallPrecisionAndF1(dictStatsCompCyntenator, arguments, benchTools.Mylimits)


# print >> sys.stderr, dictStatsCompCyntenator.keys()
# print >> sys.stderr, dictStatsCompCyntenator[(-2, -3, 1)]['familyContent']
# print >> sys.stderr, dictStatsCompCyntenator[(-1000000000, -3, 1)]['familyContent']
# print >> sys.stderr, dictStatsCompCyntenator[(-2, -3, 2)]['familyContent']
# print >> sys.stderr, dictStatsCompCyntenator[(-1000000000, -3, 2)]['familyContent']