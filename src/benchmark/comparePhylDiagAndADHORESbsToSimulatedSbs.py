#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys

import numpy

import enum
import pickle
import getpass # https://stackoverflow.com/questions/842059/is-there-a-portable-way-to-get-the-current-username-in-python
import itertools
import collections

from utils import myDiags
from utils import myTools
from utils import myADHoRe
from utils import myPhylTree
from utils import myIntervals
from utils import myLightGenomes

# TODO When MagSimus is published uncomment the next line, remove benchTools.computePairwiseSyntenyBlocksOfMagSimus,
# remove its sub functions and update links properly.
# from libs import myBreakpointsAnalyser
import benchmarkTools as benchTools

LOAD_PRECOMPUTED_BENCHMARK = True
# 5 lines to change to adapt this script to a local installation
BIN_ADHORE = myADHoRe.PATH_ADHORE_BIN

# DATA_ADHORE = '/home/jlucas/Libs/PhylDiag/data/benchmark/adhore'
# You might want to change this folder depending on where you installed PhylDiag
# /home/<user>/Libs/PhylDiag/res/benchmark/adhore
RES_ADHORE = '/home/' + getpass.getuser() + '/Libs/PhylDiag/res/benchmark/adhore'
DATA_ADHORE = RES_ADHORE
serialisationFile = RES_ADHORE + '/serialisationPhylDiagAndADHoReBenchmarkValues'


for directory in [RES_ADHORE, DATA_ADHORE]:
    if not os.path.exists(directory):
        os.makedirs(directory)

__doc__ = """This script compares the sbs of PhylDiag and the base_clusters of i-ADHoRe 3.0 with the reference
conserved segments (cs) recorded during the simulation of MagSimus (with the breakpoint analyser feature)
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

        # consider the orientation of genes for the comparison ?
        ('oriented', bool, False),

        # parameters (may contain several values, varying)
        #('gapMaxs', str, '(0,1,2,3,4,5,10,15,20,50)'),
        ('gapMaxs', str, '(0,1,2,3,4,5,6,7,8,9,10,15,20,50)'),
        #('gapMaxs', str, '(5,10)'),
        #('gapMaxs', str, '(1,)'),
        #('tandemGapMaxs', tuple, (0, 10, 100)),
        #('gapMax_DiagsHp_SbsHp', int, 0),
        #('identifyMonoGenicCs', bool, False),
        ('pThresholdPhylDiag', str, 'None'),
        ('pThresholdADHORE', float, 0.001),
        ('distanceMetric', str, 'CD'),
        ('distinguishMonoGenicDiags', bool, True),
        #('identifyMicroRearr', bool, False),

        ('preComputePairwiseSbs', bool, True),

        # edit sbs found
        ('removeSbsOfLen1', bool, False),
        ('reduceToSameGeneNames', bool, False),

        #('analyseMistakesOneByOne', bool, False),
        #('mistake', str, 'Fn'),

        ('outFigureName', str, '%s/PhylDiag_ADHoRe.svg' % RES_ADHORE)
    ],
    __doc__,
    # load Only Default Options
    loadOnlyDefaultOptions=False)

# different options of PhylDiag
# variableArg = (('imr', ('imcs', 'gmmi'), 'tm'), [(False, (False, None), None)])
variableArg = (('imr', ('imcs', 'gmmi'), 'tm'), [(False, (False, None), None),
                                                 (True, (False, None), None),
                                                 (False, (True, 0), None),
                                                 (True, (True, 0), None),
                                                 (True, (True, 1), None),
                                                 (True, (True, 1), 4)])

if not LOAD_PRECOMPUTED_BENCHMARK :
    # if arguments['gapMax_DiagsHp_SbsHp'] == 'None':
    #     arguments['gapMax_DiagsHp_SbsHp'] = None
    # else:
    #     arguments['gapMax_DiagsHp_SbsHp'] = int(arguments['gapMax_DiagsHp_SbsHp'])
    if arguments['gapMaxs'] == 'None':
        arguments['gapMaxs'] = None
    else:
        arguments['gapMaxs'] = tuple(eval(arguments['gapMaxs']))
        filterType = myDiags.FilterType.__getattr__(arguments["filterType"])
    pThresholdPhylDiag = arguments['pThresholdPhylDiag'] if arguments['pThresholdPhylDiag'] != 'None' else None
    # tandemGapMaxs = arguments['tandemGapMaxs']
    # identifyMonoGenicInversionS = [False, True]
    # pThresholdS = [1.0, 0.1, 0.01, 0.001, None]
    # gapMaxMicroInvS = [0, 1, 2, 3, 4, 5]
    # distanceMetricS = ['CD', 'MD', 'DPD', 'ED']

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

    (newGenome1, newGenome2, newFamilies) = benchTools.ensureNamesDifferBetweenGenome1AndGenome2(genome1, genome2,
                                                                                                 ancGenes)

    progressBar = myTools.ProgressBar(len(arguments['gapMaxs']) * len(variableArg[1]))
    progress = 0
    for varArg in variableArg[1]:
        dictStatsCompPhylDiag[varArg] = {}
        (imr, (imcs, gmmi), tm) = varArg
        for gapMax in arguments['gapMaxs']:
            suffix = 'gapMax%s_%s_%s_%s_%s' % (gapMax, imr, imcs, gmmi, tm)
            dictStatsCompPhylDiag[varArg][gapMax] = {}
            print >> sys.stderr, "gapMax = %s" % gapMax
            progress += 1
            progressBar.printProgressIn(sys.stderr, progress)
            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(newGenome1, newGenome2, newFamilies,
                                                                filterType=filterType,
                                                                gapMax=gapMax,
                                                                distanceMetric=arguments['distanceMetric'],
                                                                tandemGapMax=arguments['tandemGapMax'],
                                                                gapMax_Diag_Sbs=gmmi,  #arguments['gapMax_Diag_Sbs'],
                                                                distinguishMonoGenicDiags=arguments['distinguishMonoGenicDiags'],
                                                                #identifyMonoGenicCs=varArg,  #arguments['identifyMonoGenicCs'],
                                                                identifyMonoGenicCs=imcs,  #arguments['identifyMonoGenicCs'],
                                                                #identifyMicroRearr=arguments['identifyMicroRearr'],
                                                                truncationMax=tm,
                                                                identifyMicroRearr=imr,
                                                                pThreshold=pThresholdPhylDiag,
                                                                verbose=True)
            coverage = myDiags.computeMeanCoverageInGenes(newGenome1, newGenome2, sbsInPairComp, families=newFamilies)
            print >> sys.stderr, 'mean coverage in genes = %.0f%%' % (coverage * 100)

            #myDiags.printSbsFile(sbsInPairComp, newGenome1, newGenome2, sys.stderr)
            phylDiagSbsGenome = myDiags.buildGenomeFromSbs(sbsInPairComp, sbLengthThreshold=None)
            (sbsTrue1, phylDiagSbsGenome) = benchTools.editSbs(sbsTrue, phylDiagSbsGenome,
                                                               removeSbsOfLen1=arguments['removeSbsOfLen1'],
                                                               reduceToSameGeneNames=arguments['reduceToSameGeneNames'])
            nbSbsPhylDiag = len(phylDiagSbsGenome.keys())
            nbFamiliesInSbsPhylDiag = len(set(phylDiagSbsGenome.getGeneNames(asA=set, checkNoDuplicates=False)))

            # analyse ancestral gene content in sbs compared to ref sbs
            if tm is None:
                print >> sys.stderr, "analyse ancestral gene content of phyldiag sbs"
                benchTools.analyseFamiliesInSyntenyBlocks(phylDiagSbsGenome, sbsTrue1)
            else:
                countFamilyInstancesInSbsPhylDiag = collections.Counter(phylDiagSbsGenome.getGeneNames(asA=list))
                assert phylDiagSbsGenome.getGeneNames(checkNoDuplicates=True), \
                    'While ancGenes should be in only one copy in sbs (because of solving overlaps (tm not None), list of 10 ancGenes in multiple copies in PhylDiag Sbs %s' % \
                    sorted([(gn, count) for (gn, count) in countFamilyInstancesInSbsPhylDiag.iteritems()],
                           key=lambda x: x[1], reverse=True)[:10]

            comp_adjacency = myIntervals.compareGenomes(phylDiagSbsGenome, sbsTrue1,
                                                        mode='adjacency',
                                                        oriented=arguments['oriented'], verbose=True)
            comp_chromExtremity = myIntervals.compareGenomes(phylDiagSbsGenome, sbsTrue1,
                                                             mode='chromExtremity',
                                                             oriented=arguments['oriented'], verbose=True)
            comp_familyContent = myIntervals.compareGenomes(phylDiagSbsGenome, sbsTrue1,
                                                          mode='geneName',
                                                          oriented=arguments['oriented'], verbose=True)
            # sTpNs = set(gn for (gn, _) in sTp) if arguments['oriented'] else sTp
            dictStatsCompPhylDiag[varArg][gapMax]['adjacency'] = comp_adjacency
            dictStatsCompPhylDiag[varArg][gapMax]['chromExtremity'] = comp_chromExtremity
            dictStatsCompPhylDiag[varArg][gapMax]['familyContent'] = comp_familyContent
            dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'] = (nbSbsPhylDiag, nbSbsTrue)
            dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'] = (nbFamiliesInSbsPhylDiag, nbAncGenesInSbsTrue)
            dictStatsCompPhylDiag[varArg][gapMax]['coverage'] = coverage
            phylDiagSbsGenome.printIn(open(RES_ADHORE + '/phylDiag_sbs_%s.list.bz2' % suffix, 'w'))

    (newGenome1, newGenome2, newFamilies) = benchTools.ensureNamesDifferBetweenGenome1AndGenome2(genome1, genome2, ancGenes)

    dictStatsCompAdhore = {}
    for gapMax in arguments['gapMaxs']:
        # ADHoRe
        if gapMax == 0:
            # ADHoRe does not allow gapMax = 0
            continue
        print >> sys.stderr, "gapMax = %s" % gapMax
        dictStatsCompAdhore[gapMax] = {}
        assert gapMax != 0

        adhoreSbsGenome = myADHoRe.launchADHoRe(newGenome1, newGenome2, newFamilies,
                                                gapMax=gapMax, tandemGapMax=arguments['tandemGapMax'], pThreshold=arguments['pThresholdADHORE'], minimalLengthForSbs=3,
                                                filterType=filterType,
                                                outADHoReChromosomes="%s/Genome.%s.Chr%s.list" % (DATA_ADHORE, '%s', '%s'),
                                                outAHoReFamilies="%s/families.csv" % DATA_ADHORE,
                                                outAHoReConfigurationFile="%s/dataset_G1_G2.ini" % DATA_ADHORE,
                                                resADHoRePath="%s" % RES_ADHORE,
                                                pathIADHOREbin=BIN_ADHORE)
        (sbsTrue2, adhoreSbsGenome) = benchTools.editSbs(sbsTrue, adhoreSbsGenome,
                                                         removeSbsOfLen1=arguments['removeSbsOfLen1'],
                                                         reduceToSameGeneNames=arguments['reduceToSameGeneNames'])

        # analyse overlaps
        nbSbsAdhore = len(adhoreSbsGenome.keys())
        nbFamiliesInSbsAdhore = len(set(adhoreSbsGenome.getGeneNames(asA=set, checkNoDuplicates=False)))
        print >> sys.stderr, "analyse ancestral genes content of adhore sbs"
        benchTools.analyseFamiliesInSyntenyBlocks(adhoreSbsGenome, sbsTrue2)

        comp_adjacency = myIntervals.compareGenomes(adhoreSbsGenome, sbsTrue2,
                                                    mode='adjacency',
                                                    oriented=arguments['oriented'], verbose=True)
        comp_chromExtremity = myIntervals.compareGenomes(adhoreSbsGenome, sbsTrue2,
                                                         mode='chromExtremity',
                                                         oriented=arguments['oriented'], verbose=True)
        comp_familyContent = myIntervals.compareGenomes(adhoreSbsGenome, sbsTrue2,
                                                      mode='geneName',
                                                      oriented=arguments['oriented'], verbose=True)
        dictStatsCompAdhore[gapMax]['adjacency'] = comp_adjacency
        dictStatsCompAdhore[gapMax]['chromExtremity'] = comp_chromExtremity
        dictStatsCompAdhore[gapMax]['familyContent'] = comp_familyContent
        dictStatsCompAdhore[gapMax]['nbSbs'] = (nbSbsAdhore, nbSbsTrue)
        dictStatsCompAdhore[gapMax]['nbFamilies'] = (nbFamiliesInSbsAdhore, nbAncGenesInSbsTrue)
        #dictStatsCompAdhore[gapMax] = (eff, (nbSbsAdhore, nbSbsTrue), (nbFamiliesInSbsAdhore, nbAncGenesInSbsTrue))
        adhoreSbsGenome.printIn(open(RES_ADHORE + '/adhore_sbs_gapMax%s.list.bz2' % gapMax, 'w)'))

    # Serialize
    pickle.dump((dictStatsCompPhylDiag, dictStatsCompAdhore), open(serialisationFile, 'w'))

else:
    # Deserialize
    (dictStatsCompPhylDiag, dictStatsCompAdhore)  = pickle.load(open(serialisationFile, 'r'))

####################################################################
# Plot Results
####################################################################
import matplotlib.pyplot as plt
def plotRecallAndPrecision(dictStatsCompPhylDiag, dictStatsCompAdhore, arguments, variableArg):
    fig, axes = plt.subplots(nrows=2, ncols=3)
    colTitles = ['extremities', 'adjacencies', 'anc. families']
    rowTitles = ['recall', 'precision', 'F1-score']
    # mapping between gapMaxs and ticksLocations
    #gapMaxs = sorted(dictStatsCompPhylDiag.values()[0])
    gapMaxs = arguments['gapMaxs']
    tickIdxByGapMaxPhylDiag = collections.OrderedDict([(gapMax, i) for i, gapMax in enumerate(gapMaxs)])
    for ax, colTitle in zip(axes[0], colTitles):
        ax.set_title(colTitle)
    for ax, rowTitle in zip(axes[:, 0], rowTitles):
        ax.set_ylabel(rowTitle, rotation=90, size='large')
    for ax in axes[0]:
        ax.set_xticks(tickIdxByGapMaxPhylDiag.values())
        ax.set_xticklabels([str(gapMax) for gapMax in tickIdxByGapMaxPhylDiag.keys()])
    for ax in axes[1]:
        ax.set_xticks(tickIdxByGapMaxPhylDiag.values())
        ax.set_xticklabels([str(gapMax) for gapMax in tickIdxByGapMaxPhylDiag.keys()])
        ax.set_xlabel('gapMax')
    # for (i, j), ax in numpy.ndenumerate(axes):
    #     ax.set_ylim((0, 1)
    #     plt.setp(ax.get_xticklabels(), rotation=-45, horizontalalignment='left')
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(10)
    #     for label in ax.xaxis.get_ticklabels()[::2]:
    #         label.set_visible(False)
    #plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')
    axes[0][0].set_ylim((0, 1))
    axes[1][0].set_ylim((0, 1))
    lineWidth = 2.5

    variableArg_str = benchTools.variableArgIntoVariableArgString(variableArg)
    # lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linecycler = itertools.cycle(lines)
    colorcycler = itertools.cycle(colors)

    modeOfComparisons = ['chromExtremity', 'adjacency', 'familyContent']
    for i, modeOfComparison in enumerate(modeOfComparisons):
        axis_recall = axes[0][i]
        axis_precision = axes[1][i]

        algo = 'PhylDiag'
        linesPhylDiag = [None] * len(variableArg[1])
        for ivar, varArg in enumerate(variableArg[1]):
            label = variableArg_str[1][ivar]
            label = algo + ' (' + ','.join(label) + ')'
            color = next(colorcycler)
            lineStyle = next(linecycler)
            recalls = [dictStatsCompPhylDiag[varArg][gapMax][modeOfComparison].r for gapMax in gapMaxs]
            precisions = [dictStatsCompPhylDiag[varArg][gapMax][modeOfComparison].p for gapMax in gapMaxs]
            Yss = (recalls, precisions)
            for ii, (ax, Ys) in enumerate(zip((axis_recall, axis_precision), Yss)):
                #minY[subTitle] = min(Ys) if min(Ys) < minY[subTitle] else minY[subTitle]
                linesPhylDiag[ivar], = ax.plot(tickIdxByGapMaxPhylDiag.values(), Ys, color=color, linestyle=lineStyle, label=label, linewidth=lineWidth)

            # nbSbsPhylDiags = [dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'] for gapMax in gapMaxs]
            # nbFamiliesInSbsPhylDiag = [dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'] for gapMax in gapMaxs]

        # nbSbsTrues = [dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'][1] for gapMax in gapMaxs]
        # nbAncGenesInSbsTrues = [dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'][1] for gapMax in gapMaxs]

        algo = 'ADHoRe'
        label = algo
        color = next(colorcycler)
        lineStyle = next(linecycler)
        gapMaxsAdhore = sorted([gapMax for gapMax in dictStatsCompAdhore.keys() if gapMax in arguments['gapMaxs']])
        gapMaxsAdhore.remove(0)
        recalls = [dictStatsCompAdhore[gapMax][modeOfComparison].r for gapMax in gapMaxsAdhore]
        precisions = [dictStatsCompAdhore[gapMax][modeOfComparison].p for gapMax in gapMaxsAdhore]
        Yss = (recalls, precisions)
        # no gapMax = 0 for ADHoRe
        tickIdxByGapMaxADHoRe = collections.OrderedDict([(gapMax, tickIdxByGapMaxPhylDiag[gapMax]) for gapMax in gapMaxsAdhore])
        for ax, Ys in zip((axis_recall, axis_precision), Yss):
            #minY[subTitle] = min(Ys) if min(Ys) < minY[subTitle] else minY[subTitle]
            lineADhoRe, = ax.plot(tickIdxByGapMaxADHoRe.values(), Ys, color=color, linestyle=lineStyle, label=label, linewidth=lineWidth)

        # nbSbsAdhores = [dictStatsCompAdhore[gapMax]['nbSbs'][1] for gapMax in gapMaxsAdhore]
        # nbFamiliesInSbsAdhores = [dictStatsCompAdhore[gapMax]['nbFamilies'][0] for gapMax in gapMaxsAdhore]


    # change the format of y values
    for (i, j), ax in numpy.ndenumerate(axes):
        vals = ax.get_yticks()
        if not (i == 1 and j == 2):
            ax.set_yticklabels(['{:3.1f}%'.format(x * 100) for x in vals])
        else:
            ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    # Legend
    lines = linesPhylDiag + [lineADhoRe] # + [lineTrue]
    labelsPhylDiag = ['PhylDiag (' + ','.join(vs) + ')' for vs in variableArg_str[1]]
    labels = labelsPhylDiag + ['i-ADHoRe 3.0'] #+ ['Simulation']
    titleLegend = '(' + ','.join(variableArg_str[0]) + ')'
    assert len(lines) == len(labels)
    fig.legend(lines, labels, ncol=4, title=titleLegend, loc='upper center', fontsize=15)
    fig.tight_layout()
    plt.show(block=False)
    plt.savefig(arguments['outFigureName'], format='svg')
    #plt.close()
def plotF1(dictStatsCompPhylDiag, dictStatsCompAdhore, arguments, variableArg):
    fig, axes = plt.subplots(nrows=1, ncols=3, squeeze=False)
    colTitles = ['extremities', 'adjacencies', 'anc. families']
    rowTitles = ['F1-score']
    # mapping between gapMaxs and ticksLocations
    #gapMaxs = sorted(dictStatsCompPhylDiag.values()[0])
    gapMaxs = arguments['gapMaxs']
    tickIdxByGapMaxPhylDiag = collections.OrderedDict([(gapMax, i) for i, gapMax in enumerate(gapMaxs)])
    for ax, colTitle in zip(list(axes[0]), colTitles):
        ax.set_title(colTitle)
    for ax, rowTitle in zip(axes[:, 0], rowTitles):
        ax.set_ylabel(rowTitle, rotation=90, size='large')
    for ax in axes[0]:
        ax.set_xticks(tickIdxByGapMaxPhylDiag.values())
        ax.set_xticklabels([str(gapMax) for gapMax in tickIdxByGapMaxPhylDiag.keys()])
    for ax in axes[0]:
        ax.set_xticks(tickIdxByGapMaxPhylDiag.values())
        ax.set_xticklabels([str(gapMax) for gapMax in tickIdxByGapMaxPhylDiag.keys()])
        ax.set_xlabel('gapMax')
    # for (i, j), ax in numpy.ndenumerate(axes):
    #     ax.set_ylim((0, 1)
    #     plt.setp(ax.get_xticklabels(), rotation=-45, horizontalalignment='left')
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(10)
    #     for label in ax.xaxis.get_ticklabels()[::2]:
    #         label.set_visible(False)
    #plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')
    axes[0][0].set_ylim((0, 1))
    #axes[1][0].set_ylim((0, 1))
    lineWidth = 2.5

    variableArg_str = benchTools.variableArgIntoVariableArgString(variableArg)
    # lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linecycler = itertools.cycle(lines)
    colorcycler = itertools.cycle(colors)

    modeOfComparisons = ['chromExtremity', 'adjacency', 'familyContent']
    for i, modeOfComparison in enumerate(modeOfComparisons):
        axis_F1 = axes[0][i]

        algo = 'PhylDiag'
        linesPhylDiag = [None] * len(variableArg[1])
        for ivar, varArg in enumerate(variableArg[1]):
            label = variableArg_str[1][ivar]
            label = algo + ' (' + ','.join(label) + ')'
            color = next(colorcycler)
            lineStyle = next(linecycler)
            F1s = [dictStatsCompPhylDiag[varArg][gapMax][modeOfComparison].f1 for gapMax in gapMaxs]
            linesPhylDiag[ivar], = axis_F1.plot(tickIdxByGapMaxPhylDiag.values(), F1s, color=color, linestyle=lineStyle, label=label, linewidth=lineWidth)

            # nbSbsPhylDiags = [dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'] for gapMax in gapMaxs]
            # nbFamiliesInSbsPhylDiag = [dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'] for gapMax in gapMaxs]

        # nbSbsTrues = [dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'][1] for gapMax in gapMaxs]
        # nbAncGenesInSbsTrues = [dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'][1] for gapMax in gapMaxs]

        algo = 'ADHoRe'
        label = algo
        color = next(colorcycler)
        lineStyle = next(linecycler)
        gapMaxsAdhore = sorted([gapMax for gapMax in dictStatsCompAdhore.keys() if gapMax in arguments['gapMaxs']])
        gapMaxsAdhore.remove(0)
        F1s = [dictStatsCompAdhore[gapMax][modeOfComparison].f1 for gapMax in gapMaxsAdhore]
        # no gapMax = 0 for ADHoRe
        tickIdxByGapMaxADHoRe = collections.OrderedDict([(gapMax, tickIdxByGapMaxPhylDiag[gapMax]) for gapMax in gapMaxsAdhore])
        lineADhoRe, = axis_F1.plot(tickIdxByGapMaxADHoRe.values(), F1s, color=color, linestyle=lineStyle, label=label, linewidth=lineWidth)

        # nbSbsAdhores = [dictStatsCompAdhore[gapMax]['nbSbs'][1] for gapMax in gapMaxsAdhore]
        # nbFamiliesInSbsAdhores = [dictStatsCompAdhore[gapMax]['nbFamilies'][0] for gapMax in gapMaxsAdhore]


    # change the format of y values
    for (i, j), ax in numpy.ndenumerate(axes):
        vals = ax.get_yticks()
        if not (i == 1 and j == 2):
            ax.set_yticklabels(['{:3.1f}%'.format(x * 100) for x in vals])
        else:
            ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])

    # Legend
    lines = linesPhylDiag + [lineADhoRe] # + [lineTrue]
    labelsPhylDiag = ['PhylDiag (' + ','.join(vs) + ')' for vs in variableArg_str[1]]
    labels = labelsPhylDiag + ['i-ADHoRe 3.0'] #+ ['Simulation']
    titleLegend = '(' + ','.join(variableArg_str[0]) + ')'
    assert len(lines) == len(labels)
    fig.legend(lines, labels, ncol=4, title=titleLegend, loc='upper center', fontsize=15)
    fig.tight_layout()
    plt.show(block=True)
    plt.savefig(arguments['outFigureName'], format='svg')
    plt.close()

def plotRecallPrecisionAndF1(dictStatsCompPhylDiag, dictStatsCompAdhore, arguments, variableArg, Mylimits):
    fig, axes = plt.subplots(nrows=3, ncols=3)
    colTitles = ['extremities', 'adjacencies', 'gene names']
    rowTitles = ['recall', 'precision', 'F1-score']
    # mapping between gapMaxs and ticksLocations
    #gapMaxs = sorted(dictStatsCompPhylDiag.values()[0])
    gapMaxs = arguments['gapMaxs']
    tickIdxByGapMaxPhylDiag = collections.OrderedDict([(gapMax, i) for i, gapMax in enumerate(gapMaxs)])
    for ax, colTitle in zip(axes[0], colTitles):
        ax.set_title(colTitle)
    for ax, rowTitle in zip(axes[:, 0], rowTitles):
        ax.set_ylabel(rowTitle, rotation=90, size='large')
    for ax in axes.flatten():
        ax.set_xticks(tickIdxByGapMaxPhylDiag.values())
        ax.set_xticklabels([str(gapMax) for gapMax in tickIdxByGapMaxPhylDiag.keys()])
    for ax in axes[2]:
        ax.set_xlabel('gapMax')
    # for (i, j), ax in numpy.ndenumerate(axes):
    #     ax.set_ylim((0, 1)
    #     plt.setp(ax.get_xticklabels(), rotation=-45, horizontalalignment='left')
    #     for tick in ax.xaxis.get_major_ticks():
    #         tick.label.set_fontsize(10)
    #     for label in ax.xaxis.get_ticklabels()[::2]:
    #         label.set_visible(False)
    #plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

    # axes limits
    # axes[0][0].set_ylim((0, 1))
    # axes[1][0].set_ylim((0, 1))
    # axes[2][0].set_ylim((0, 1))
    for (i, j), yl in numpy.ndenumerate(Mylimits):
        axes[i][j].set_ylim(yl)
    lineWidth = 2.5

    variableArg_str = benchTools.variableArgIntoVariableArgString(variableArg)
    # lines = ["-", "--", "-.", ":"]
    lines = ["-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linecycler = itertools.cycle(lines)
    colorcycler = itertools.cycle(colors)

    modeOfComparisons = ['chromExtremity', 'adjacency', 'familyContent']
    for i, modeOfComparison in enumerate(modeOfComparisons):
        axis_recall = axes[0][i]
        axis_precision = axes[1][i]
        axis_f1 = axes[2][i]

        algo = 'PhylDiag'
        linesPhylDiag = [None] * len(variableArg[1])
        for ivar, varArg in enumerate(variableArg[1]):
            label = variableArg_str[1][ivar]
            label = algo + ' (' + ','.join(label) + ')'
            color = next(colorcycler)
            lineStyle = next(linecycler)
            recalls = [dictStatsCompPhylDiag[varArg][gapMax][modeOfComparison].r for gapMax in gapMaxs]
            precisions = [dictStatsCompPhylDiag[varArg][gapMax][modeOfComparison].p for gapMax in gapMaxs]
            f1s = [dictStatsCompPhylDiag[varArg][gapMax][modeOfComparison].f1 for gapMax in gapMaxs]
            Yss = (recalls, precisions, f1s)
            for ii, (ax, Ys) in enumerate(zip((axis_recall, axis_precision, axis_f1), Yss)):
                #minY[subTitle] = min(Ys) if min(Ys) < minY[subTitle] else minY[subTitle]
                linesPhylDiag[ivar], = ax.plot(tickIdxByGapMaxPhylDiag.values(), Ys, color=color, linestyle=lineStyle, label=label, linewidth=lineWidth)

            # nbSbsPhylDiags = [dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'] for gapMax in gapMaxs]
            # nbFamiliesInSbsPhylDiag = [dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'] for gapMax in gapMaxs]

        # nbSbsTrues = [dictStatsCompPhylDiag[varArg][gapMax]['nbSbs'][1] for gapMax in gapMaxs]
        # nbAncGenesInSbsTrues = [dictStatsCompPhylDiag[varArg][gapMax]['nbFamilies'][1] for gapMax in gapMaxs]

        algo = 'ADHoRe'
        label = algo
        color = next(colorcycler)
        lineStyle = next(linecycler)
        gapMaxsAdhore = sorted([gapMax for gapMax in dictStatsCompAdhore.keys() if gapMax in arguments['gapMaxs']])
        if 0 in gapMaxsAdhore:
            gapMaxsAdhore.remove(0)
        recalls = [dictStatsCompAdhore[gapMax][modeOfComparison].r for gapMax in gapMaxsAdhore]
        precisions = [dictStatsCompAdhore[gapMax][modeOfComparison].p for gapMax in gapMaxsAdhore]
        f1s = [dictStatsCompAdhore[gapMax][modeOfComparison].f1 for gapMax in gapMaxsAdhore]
        Yss = (recalls, precisions, f1s)
        # no gapMax = 0 for ADHoRe
        tickIdxByGapMaxADHoRe = collections.OrderedDict([(gapMax, tickIdxByGapMaxPhylDiag[gapMax]) for gapMax in gapMaxsAdhore])
        for ax, Ys in zip((axis_recall, axis_precision, axis_f1), Yss):
            #minY[subTitle] = min(Ys) if min(Ys) < minY[subTitle] else minY[subTitle]
            lineADhoRe, = ax.plot(tickIdxByGapMaxADHoRe.values(), Ys, color=color, linestyle=lineStyle, label=label, linewidth=lineWidth)

        # nbSbsAdhores = [dictStatsCompAdhore[gapMax]['nbSbs'][1] for gapMax in gapMaxsAdhore]
        # nbFamiliesInSbsAdhores = [dictStatsCompAdhore[gapMax]['nbFamilies'][0] for gapMax in gapMaxsAdhore]


    # change the format of x and y values
    #setOfAxesCoordsWithIncreasedPrecision = {(1,2), (2,2)}
    #setOfAxesCoordsWithIncreasedPrecision = {(1,2)}
    setOfAxesCoordsWithIncreasedPrecision = {}
    for (i, j), ax in numpy.ndenumerate(axes):
        ax.tick_params(labelsize=10)
        vals = ax.get_yticks()
        if not (i,j) in setOfAxesCoordsWithIncreasedPrecision:
            ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in vals])
        else:
            ax.set_yticklabels(['{:3.1f}%'.format(x * 100) for x in vals])

    # Legend
    lines = linesPhylDiag + [lineADhoRe] # + [lineTrue]
    labelsPhylDiag = ['PhylDiag (' + ','.join(vs) + ')' for vs in variableArg_str[1]]
    labels = labelsPhylDiag + ['i-ADHoRe 3.0'] #+ ['Simulation']
    titleLegend = '(' + ','.join(variableArg_str[0]) + ')'
    assert len(lines) == len(labels)
    fig.legend(lines, labels, ncol=4, title=titleLegend, loc='upper center', fontsize=12)
    fig.tight_layout()
    plt.show(block=True)
    plt.savefig(arguments['outFigureName'], format='svg')
    #plt.close()

arguments['gapMaxs'] = (0,1,2,3,4,5,10,15,20,50)
plotRecallPrecisionAndF1(dictStatsCompPhylDiag, dictStatsCompAdhore, arguments, variableArg, benchTools.Mylimits)