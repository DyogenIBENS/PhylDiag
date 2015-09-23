#!/usr/bin/python
# -*- coding: utf-8 -*-

#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France
from utils import myMapping

__doc__ = """
    Plot the distribution of the lengths of synteny blocks (either in genes or in bp)
"""

import sys

from matplotlib import pyplot as plt
import numpy as np

import utils.myTools as myTools
import utils.myDiags as myDiags
import utils.myGenomes as myGenomes
import utils.myLightGenomes as myLightGenomes
import scipy.optimize
import random
from scipy.stats import expon

def Exp(t, theta):
    return 1.0/theta * np.exp(-t/theta)

def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

import scipy as sp
def fit_exp_nonlinear(t, y, p0=(170, -0.1, 0)):
    (opt_parms, parm_cov) = sp.optimize.curve_fit(model_func, t, y, p0=p0)
    A, K, C = opt_parms
    return A, K, C

def fit_exp_linear(t, y, C=0, removeNullValues=True, maxT=None, doPlot=True):
    assert len(t) == len(y)
    if removeNullValues:
        # remove all values <= 0
        newt = []
        newy = []
        for (vt, vy) in zip(t, y):
            if vy > 0:
                if maxT is None or vt <= maxT:
                    newt.append(vt)
                    newy.append(vy)
    t = np.asarray(newt)
    y = np.asarray(newy)
    assert len(t) == len(y)
    y = y - C
    # y = A * exp(Kx)
    y = np.log(y)
    # log(y) = A_log + K*x
    K, A_log = np.polyfit(x=t, y=y, deg=1)
    if doPlot:
        plt.figure()
        plt.plot(t, y, marker='o', linestyle='', color='k')
        plt.plot(t, K * t + A_log, marker='o', linestyle='', color='b')
        plt.xlabel('lengths of sbs (logarithmic scale)')
        plt.ylabel('#sbs')
        # plt.show()
    A = np.exp(A_log)
    return A, K

def lengthProjectionOnGenome(rankGenome, sb, c, genome, correctLens=True, meanInterGeneLen=0.0):
    lGeneCoord = []
    sblX = sb.l1 if rankGenome == 1 else sb.l2
    for tbX in sblX:
        for gIdx in tbX:
            #g1Pos = genome1.getPosition([g1n]).pop()
            #chromosome = g1Pos.chromosome
            #assert c1 == chromosome
            #index = g1Pos.index
            c = myGenomes.commonChrName(c)
            g = genome.lstGenes[c][gIdx]
            lGeneCoord.extend([g.beginning, g.end])
    minOnG = min(lGeneCoord)
    maxOnG = max(lGeneCoord)
    lengthG = maxOnG - minOnG
    if correctLens:
        # see Nadeau & Taylor 1984 815-816
        n = len(sblX)
        r = lengthG
        if n >= 2:
            # this consider the local density of markers
            m = r * float(n+1)/float(n-1)
        elif n == 1:
            # add for each extremity (twice), the mean intergene length / 2
            m = r + 2 * float(meanInterGeneLen) / 2
        else:
            raise ValueError('A sb should at least contain one marker')
        lengthG = m
    return lengthG

def plotDensityFunction(ax, lSbsLengths, bins, binwidth, minSbLength):
    print >> sys.stderr, '\n'.join([str(lSbsLengths), str(bins), str(binwidth), str(minSbLength)])
    n, _, patches = ax.hist(lSbsLengths, bins=bins, histtype='bar')
    t = np.asarray(bins[1:])
    v = np.asarray(n)
    print >> sys.stderr, 'meanSbLength =', float(sum(lSbsLengths)) / len(lSbsLengths)

    def computeMeanSbLength(lSbsLengths, binwidth=None, minSbLength=None):
        # first fit an exponential of the form s * 1/theta * exp(-x/theta)
        # with s the scale factor = total #sbs
        loc, scale = expon.fit(lSbsLengths, floc=0.0)
        assert loc == 0.0

        #A, K = fit_exp_linear(t, v, 0)  # A exp(K t)
        # A = s * 1/theta
        #theta = -1.0/K

        # proportion of the missing normalised distribution between 0 and minSbLength
        # propMissing = 1.0 / np.exp(-minSbLength/scale) - 1.0
        # s = (propMissing + 1) * len(lSbsLengths)
        int_0_removedHeadLength_expon = np.exp(- float(minSbLength) / scale)
        s = len(lSbsLengths) / float(int_0_removedHeadLength_expon)

        theta = scale
        # print >> sys.stderr, "theta=", theta
        # print >> sys.stderr, "s=", s
        missingBins = list(np.arange(0, int(minSbLength), binwidth))
        nbMissingSbs = s * Exp(np.asarray([bin - binwidth/2.0 for bin in missingBins]), theta)
        missingNbSbsPerBin = zip(missingBins, nbMissingSbs)
        lSbsLengthsMissing = [(bin + random.random() * binwidth) for (bin, sbLength) in missingNbSbsPerBin for _ in range(int(sbLength))]
        lSbsLengthsc = list(lSbsLengthsMissing) + lSbsLengths
        meanSbLength = float(sum(lSbsLengthsc) / len(lSbsLengthsc))
        # TODO compare s and len(lSbsLengthsc)
        print >> sys.stderr, "meanSbLength=%s, scale=%s" % (meanSbLength, theta)
        print >> sys.stderr, "s=%s, len(lSbsLengthsc)=%s" % (s, len(lSbsLengthsc))
        return (meanSbLength, lSbsLengthsc, missingNbSbsPerBin)

    (meanSbLengthc, lSbsLengthsc, missingNbSbsPerBin) = computeMeanSbLength(lSbsLengths, binwidth=binwidth, minSbLength=minSbLength)
    print >> sys.stderr, 'meanSbLength (corrected) =', meanSbLengthc
    # negatie Exponential Distribution Parameterisation
    theta = meanSbLength
    ysModel = Exp(t, theta)
    # from normal distribuion (area = 1) to a distribution of area = nb of sbs
    ysModel = [y*len(lSbsLengths) for y in ysModel]
    thetac = meanSbLengthc
    ysModelc = Exp(t, thetac)
    # from normal distribuion (area = 1) to a distribution of area = nb of sbs
    ysModelc = [y*len(lSbsLengthsc) for y in ysModelc]

    ax.plot(bins[1:], ysModel, color='r', linewidth=1.0, label="$y = %s \\times \\frac{1}{%0.2f} e^{\\frac{-t}{%0.2f}}$" % (len(lSbsLengths), theta, theta), zorder=100)
    ax.plot(bins[1:], ysModelc, color='k', linewidth=1.0, label="with corrected $\\theta$, $y = %s \\times \\frac{1}{%0.2f} e^{\\frac{-t}{%0.2f}}$" % (len(lSbsLengthsc), thetac, thetac), zorder=100)
    for (bin, nbSbsMissing) in missingNbSbsPerBin:
        ax.bar(bin, nbSbsMissing, binwidth, color='g', zorder=0)
    ax.set_ylabel('Nb of synteny blocks')
    # str_maxSbLength = "\infty" if maxSbLength > 100 else maxSbLength
    # ax.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], minSbLength, str_maxSbLength))
    ax.legend()

def plotDensityFunction2(ax, lSbsLengths, bins, binwidth, minSbLength):
    print >> sys.stderr, '\n'.join([str(lSbsLengths), str(bins), str(binwidth), str(minSbLength)])
    n, _, patches = ax.hist(lSbsLengths, bins=bins, histtype='bar')
    t = np.asarray(bins[1:])
    v = np.asarray(n)
    print >> sys.stderr, 'meanSbLength =', float(sum(lSbsLengths)) / len(lSbsLengths)

    def computeMeanSbLength(lSbsLengths, binwidth=None, minSbLength=None):
        # first fit an exponential of the form s * 1/theta * exp(-x/theta)
        # with s the scale factor = total #sbs
        loc, scale = expon.fit(lSbsLengths, floc=0.0)
        assert loc == 0.0

        #A, K = fit_exp_linear(t, v, 0)  # A exp(K t)
        # A = s * 1/theta
        #theta = -1.0/K

        # proportion of the missing normalised distribution between 0 and minSbLength
        # propMissing = 1.0 / np.exp(-minSbLength/scale) - 1.0
        # s = (propMissing + 1) * len(lSbsLengths)
        int_0_removedHeadLength_expon = np.exp(- float(minSbLength) / scale)
        s = len(lSbsLengths) / float(int_0_removedHeadLength_expon)

        theta = scale
        missingBins = list(np.arange(0, int(minSbLength), binwidth))
        nbMissingSbs = s * Exp(np.asarray([bin - binwidth/2.0 for bin in missingBins]), theta)
        missingNbSbsPerBin = zip(missingBins, nbMissingSbs)
        lSbsLengthsMissing = [(bin + random.random() * binwidth) for (bin, sbLength) in missingNbSbsPerBin for _ in range(int(sbLength))]
        lSbsLengthsc = list(lSbsLengthsMissing) + lSbsLengths
        meanSbLength = float(sum(lSbsLengthsc) / len(lSbsLengthsc))
        print >> sys.stderr, "meanSbLength=%s, scale=%s" % (meanSbLength, theta)
        print >> sys.stderr, "s=%s, len(lSbsLengthsc)=%s" % (s, len(lSbsLengthsc))
        return (meanSbLength, lSbsLengthsc, missingNbSbsPerBin)

    (meanSbLengthc, lSbsLengthsc, missingNbSbsPerBin) = computeMeanSbLength(lSbsLengths, binwidth=binwidth, minSbLength=minSbLength)
    print >> sys.stderr, 'meanSbLength (corrected) =', meanSbLengthc
    # negatie Exponential Distribution Parameterisation
    theta = meanSbLength
    ysModel = Exp(t, theta)
    # from normal distribuion (area = 1) to a distribution of area = nb of sbs
    ysModel = [y*len(lSbsLengths) for y in ysModel]
    thetac = meanSbLengthc
    ysModelc = Exp(t, thetac)
    # from normal distribuion (area = 1) to a distribution of area = nb of sbs
    ysModelc = [y*len(lSbsLengthsc) for y in ysModelc]

    ax.plot(bins[1:], ysModel, color='r', linewidth=1.0, label="$y = %s \\times \\frac{1}{%0.2f} e^{\\frac{-t}{%0.2f}}$" % (len(lSbsLengths), theta, theta), zorder=100)
    ax.plot(bins[1:], ysModelc, color='k', linewidth=1.0, label="with corrected $\\theta$, $y = %s \\times \\frac{1}{%0.2f} e^{\\frac{-t}{%0.2f}}$" % (len(lSbsLengthsc), thetac, thetac), zorder=100)
    for (bin, nbSbsMissing) in missingNbSbsPerBin:
        ax.bar(bin, nbSbsMissing, binwidth, color='g', zorder=0)
    ax.set_ylabel('Nb of synteny blocks')
    # str_maxSbLength = "\infty" if maxSbLength > 100 else maxSbLength
    # ax.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], minSbLength, str_maxSbLength))
    ax.legend()

def computeMeanInterHomologLen(genome, setConsideredGeneNames):
    assert isinstance(genome, myGenomes.Genome)
    chromosomesMeanInterConsideredGeneLen = {}
    lGeneCoord = []
    totalLenChrom = 0
    totalLenConsideredGenes = 0
    for c in genome.chrList[myGenomes.ContigType.Chromosome]:
        lenConsideredGenes = 0
        nbConsideredGenes = 0
        for g in genome.lstGenes[c]:
            assert g.beginning < g.end
            lGeneCoord.extend([g.beginning, g.end])
            assert len(g.names) == 1, g.names
            if g.names[0] in setConsideredGeneNames:
                lenConsideredGenes += g.end - g.beginning
                nbConsideredGenes += 1
        # convert in Mb
        minC = float(min(lGeneCoord)) / 1000000
        maxC = float(max(lGeneCoord)) / 1000000
        lenChrom = maxC - minC
        print >> sys.stderr, c
        chromosomesMeanInterConsideredGeneLen[str(c)] = float(lenChrom) / lenConsideredGenes if lenConsideredGenes > 0 else float(lenChrom)
        totalLenChrom += lenChrom if lenChrom > 0 else None
        totalLenConsideredGenes += lenConsideredGenes
    if totalLenChrom:
        densityInConsideredGenes = totalLenConsideredGenes / totalLenChrom
    else:
        densityInConsideredGenes = None
    return (chromosomesMeanInterConsideredGeneLen, densityInConsideredGenes)

# Arguments
modesOrthos = list(myDiags.FilterType._keys)
arguments = myTools.checkArgs(
    [
        ('syntenyBlocks', file),
        ('genome1', file),
        ('genome2', file)
    ],
    [
        ('lengthUnit', str, 'Mb'),
        ('minSbLength', float, 1),
        ('maxSbLength', str, 'None')  # this should not be changed, otherwise the computation of the meanInterval is wrong
    ],
    __doc__)
minSbLength = arguments['minSbLength']
if arguments['maxSbLength'] == 'None':
    maxSbLength = sys.maxint
else:
    try:
        maxSbLength = int(arguments['maxSbLength'])
    except:
        raise ValueError('maxSbLength should is not an integer')
lengthUnit = arguments['lengthUnit']
assert lengthUnit == 'Mb' or lengthUnit == 'gene'

# load data
genome1 = myGenomes.Genome(arguments["genome1"])
genome1L = myLightGenomes.LightGenome(genome1)
print >> sys.stderr, "Genome1"
print >> sys.stderr, "nb of Chr = ", len(genome1L)
genome2 = myGenomes.Genome(arguments["genome2"])
genome2L = myLightGenomes.LightGenome(genome2)
print >> sys.stderr, "Genome2"
print >> sys.stderr, "nb of Chr = ", len(genome2L)
sbsInPairComp = myDiags.parseSbsFile(arguments['syntenyBlocks'], genome1=genome1L, genome2=genome2L)

# preprocess data
# record all mapped homologs
setOfMappedHomologs1 = set()
setOfMappedHomologs2 = set()
for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
    sblX = sb.l1
    for tbX in sblX:
        for gIdx in tbX:
            g = genome1L[c1][gIdx]
            setOfMappedHomologs1.add(g.n)
    sblX = sb.l2
    for tbX in sblX:
        for gIdx in tbX:
            g = genome2L[c2][gIdx]
            setOfMappedHomologs2.add(g.n)
lSbsLengths = []
# compute sb lengths
if arguments['lengthUnit'] == 'gene':
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        nbAncGenes = len(sb.la)
        lSbsLengths.append(nbAncGenes)
elif arguments['lengthUnit'] == 'Mb':
    (meanInterHomologLenInChr1, densityInG1) = computeMeanInterHomologLen(genome1, setOfMappedHomologs1)
    (meanInterHomologLenInChr2, densityInG2) = computeMeanInterHomologLen(genome2, setOfMappedHomologs2)
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        lengthG1 = lengthProjectionOnGenome(1, sb, c1, genome1, correctLens=True, meanInterGeneLen=meanInterHomologLenInChr1[c1])
        lengthG2 = lengthProjectionOnGenome(2, sb, c2, genome2, correctLens=True, meanInterGeneLen=meanInterHomologLenInChr2[c2])
        averageSbLength = float(lengthG1 + lengthG2) / 2.0
        # in megabases
        averageSbLength = averageSbLength / 1000000
        lSbsLengths.append(averageSbLength)
lSbsLengths.sort()
#print >> sys.stderr, lSbsLengths

# remove lengths not in [minSbLength, maxSbLength]
lSbsLengthsInInterval = [sbLength for sbLength in lSbsLengths if minSbLength <= sbLength <= maxSbLength]
print >> sys.stderr, "nb of removed sb not in interval = %s" % (len(lSbsLengths) - len(lSbsLengthsInInterval))

##########################
# Graph Parameters
##########################
binwidth=1.0
mindata=0
maxdata=20
#maxdata=max(lSbsLengths)
##########################

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(15, 6))
bins = list(np.arange(mindata, maxdata + binwidth, binwidth))
t = np.asarray(bins[1:])
##########################################################
# graphs with all sbLength in [minSbLength, maxSbLength]
##########################################################
# Since all sbs of lengths < minSbLengths are removed the next computation is biased
# meanSbLength = float(sum(lSbsLengthsInInterval))/len(lSbsLengthsInInterval)
# The mean needs to take into account the fact that only sbs of at least minSbLengths are considered
meanSbLength = float(sum(lSbsLengthsInInterval))/len(lSbsLengthsInInterval)

# equal to 8.09 Mb ~= 8 cM "One centiMorgan corresponds to about 1 million base pairs in humans on average"
print >> sys.stderr, "meanSbLength (inInterVal) = %s %s" % (meanSbLength, lengthUnit)
plotDensityFunction(ax1, lSbsLengthsInInterval, bins, binwidth, minSbLength)
ax1.set_xlim((0, maxdata))

# cumulative distribution
print >> sys.stderr, bins[1:]
n, _, patches = ax3.hist(lSbsLengthsInInterval, bins=bins, normed=1, histtype='step', cumulative=True)
# add the theoretical cumulative function
# t = np.asarray(bins[1:])
#noisy = np.asarray(n)

ys = 1 - np.exp(-t/meanSbLength)
ax3.plot(t, ys)
ax3.set_xlim((0, maxdata))
ax3.set_ylim((0, 1))
ax3.set_ylabel("% synteny blocks")
str_maxSbLength = "\infty" if maxSbLength > 100 else maxSbLength
ax3.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], minSbLength, str_maxSbLength))
ax3.legend()

#############################
# graphs with all sbLengths
#############################
meanSbLength = float(sum(lSbsLengths))/len(lSbsLengths)
plotDensityFunction(ax2, lSbsLengths, bins, binwidth, 0)
ax2.set_xlim((0, maxdata))

# cumulative distribution
n, _, patches = ax4.hist(lSbsLengths, bins=bins, normed=1, histtype='step', cumulative=True)
# add the theoretical cumulative function
# t = np.asarray(bins[1:])
#noisy = np.asarray(n)
ys = 1 - np.exp(-t/meanSbLength)
ax4.plot(t, ys)
ax4.set_xlim((0, maxdata))
ax4.set_ylim((0, 1))
ax4.set_ylabel("% synteny blocks")
ax4.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], 0, "\infty"))
ax4.legend()

fig.suptitle("Distribution of the lengths of synteny blocks\n%s-%s comparison" % (genome1.name, genome2.name))
plt.show()

##################
# Some tests
##################
def testScipyExponential():
    data0 = expon.rvs(scale=10, size=1000)
    ###################
    data = data0
    plt.figure()
    x = np.linspace(0, 100, 100)
    plt.hist(data, bins=x, normed=True)
    plt.plot(x, expon.pdf(x, loc=0, scale=10), color='g')

    #loc, scale = expon.fit(data, floc=0)
    #plt.plot(x, expon.pdf(x, loc=loc, scale=scale), color='r')

    removedHeadLength = 1.0
    dataNoHead = [v for v in data if v > removedHeadLength]
    loc1, scale1 = expon.fit(dataNoHead)
    plt.plot(x, expon.pdf(x, loc=0, scale=scale1), color='b')
    loc, scale = expon.fit(dataNoHead, floc=removedHeadLength)
    plt.plot(x, expon.pdf(x, loc=0, scale=scale), color='r')

    # non-normed graphs
    # plt.figure()
    # plt.hist(data0, bins=x, normed=False)
    # plt.plot(x, expon.pdf(x, loc=0, scale=10)*len(data0), color='r')

    plt.figure()
    plt.hist(dataNoHead, bins=x, normed=False)
    # s = len(dataNoHead) / sInvNormalisation = intergral(pdf, removedHeadLength, infty)
    int_0_removedHeadLength_expon = np.exp(- float(removedHeadLength) / scale)
    s = len(dataNoHead) / int_0_removedHeadLength_expon
    s0 = len(data0) / 1.0
    print >> sys.stderr, s0, len(dataNoHead), s
    plt.plot(x, expon.pdf(x, loc=0, scale=scale)*s, color='r')
##############################################################
# deprecated
##############################################################

# non-linear fit
#A, K, C = fit_exp_nonlinear(t, noisy)

# linear fit with the constant set to 0
# C = 0
# A, K = fit_exp_linear(t, noisy, C)
# ysModel = model_func(t, A, K, C)

#plt.tight_layout()
#plt.xlim(0, 100)
#plt.title("OSB length distribution in Human-Mouse comparison \n Confidence Interval : %s*sigma around mean" % arguments["ICfactorOfSigma"] )
#plt.title("")
# #plt.legend()
#plt.savefig(sys.stdout, format='svg')