#!/usr/bin/python
# -*- coding: utf-8 -*-

#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

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

def Exp(t, theta):
    return 1.0/theta * np.exp(- t/theta)

def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

import scipy as sp
def fit_exp_nonlinear(t, y, p0=(170, -0.1, 0)):
    (opt_parms, parm_cov) = sp.optimize.curve_fit(model_func, t, y, p0=p0)
    A, K, C = opt_parms
    return A, K, C

def fit_exp_linear(t, y, C=0, removeNullValues=True, maxT=None):
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
    y = np.log(y)
    K, A_log = np.polyfit(x=t, y=y, deg=1)
    plt.figure()
    plt.plot(t, y, marker='o', linestyle='', color='k')
    plt.plot(t, K * t + A_log, marker='o', linestyle='', color='b')
    plt.xlabel('lengths of sbs (logarithmic scale)')
    plt.ylabel('#sbs')
    A = np.exp(A_log)
    return A, K

def lengthProjectionOnGenome(rankGenome, sb, c, genome):
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
    return lengthG

def plotDensityFunction(ax, lSbsLengths, bins):
    n, _, patches = ax.hist(lSbsLengths, bins=bins, histtype='bar')
    t = np.asarray(bins[1:])
    #noisy = np.asarray(n)

    # negatie Exponential Distribution Parameterisation
    meanSbLength = float(sum(lSbsLengths))/len(lSbsLengths)
    # equal to 8.09 Mb ~= 8 cM "One centiMorgan corresponds to about 1 million base pairs in humans on average"
    print >> sys.stderr, "meanSbLength (inInterVal) = %s %s" % (meanSbLength, lengthUnit)
    theta = meanSbLength
    ysModel = Exp(t, theta)
    # from normal distribuion (area = 1) to a distribution of area = nb of sbs
    ysModel = [y*len(lSbsLengths) for y in ysModel]

    ax.plot(bins[1:], ysModel, color='k', linewidth=1.0, label="$y = %s \\times \\frac{1}{%0.2f} e^{\\frac{-t}{%0.2f}}$" % (len(lSbsLengths), theta, theta))
    ax.set_ylabel('Nb of synteny blocks')
    # str_maxSbLength = "\infty" if maxSbLength > 100 else maxSbLength
    # ax.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], minSbLength, str_maxSbLength))
    ax.legend()


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
        ('maxSbLength', str, 'None')
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

lSbsLengths = []
# compute sb lengths
if arguments['lengthUnit'] == 'gene':
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        nbAncGenes = len(sb.la)
        lSbsLengths.append(nbAncGenes)
elif arguments['lengthUnit'] == 'Mb':
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        lengthG1 = lengthProjectionOnGenome(1, sb, c1, genome1)
        lengthG2 = lengthProjectionOnGenome(2, sb, c2, genome2)
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
#maxdata=80
maxdata=max(lSbsLengths)
##########################

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(15, 6))
bins = list(np.arange(mindata, maxdata + binwidth, binwidth))
t = np.asarray(bins[1:])
##########################################################
# graphs with all sbLength in [minSbLength, maxSbLength]
##########################################################
plotDensityFunction(ax1, lSbsLengthsInInterval, bins)
ax1.set_xlim((0, maxdata))

# cumulative distribution
n, _, patches = ax3.hist(lSbsLengthsInInterval, bins=bins, normed=1, histtype='step', cumulative=True)
# add the theoretical cumulative function
# t = np.asarray(bins[1:])
#noisy = np.asarray(n)
meanSbLength = float(sum(lSbsLengthsInInterval))/len(lSbsLengthsInInterval)
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
plotDensityFunction(ax2, lSbsLengths, bins)
ax2.set_xlim((0, maxdata))

# cumulative distribution
n, _, patches = ax4.hist(lSbsLengths, bins=bins, normed=1, histtype='step', cumulative=True)
# add the theoretical cumulative function
# t = np.asarray(bins[1:])
#noisy = np.asarray(n)
meanSbLength = float(sum(lSbsLengths))/len(lSbsLengths)
ys = 1 - np.exp(-t/meanSbLength)
ax4.plot(t, ys)
ax4.set_xlim((0, maxdata))
ax4.set_ylim((0, 1))
ax4.set_ylabel("% synteny blocks")
ax4.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], 0, "\infty"))
ax4.legend()

fig.suptitle("Distribution of the lengths of synteny blocks\nHuman-Mouse comparison")
plt.show()

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