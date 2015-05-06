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

def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

import scipy as sp
def fit_exp_nonlinear(t, y, p0=(170, -0.1, 0)):
    (opt_parms, parm_cov) = sp.optimize.curve_fit(model_func, t, y, p0=p0)
    A, K, C = opt_parms
    return A, K, C

def fit_exp_linear(t, y, C=0):
    assert len(t) == len(y)
    # remove all values <= 0
    newt = []
    newy = []
    for (vt, vy) in zip(t, y):
        if vy > 0:
            newt.append(vt)
            newy.append(vy)
    t = np.asarray(newt)
    y = np.asarray(newy)
    assert len(t) == len(y)
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(x=t, y=y, deg=1)
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
        ('minShownLength', float, 1),
        ('maxShownLength', str, 'None')
    ],
    __doc__)

minShownLength = arguments['minShownLength']
if arguments['maxShownLength'] == 'None':
    maxShownLength = sys.maxint
else:
    try:
        maxShownLength = int(arguments['maxShownLength'])
    except:
        raise ValueError('maxShownLength should is not an integer')

lengthUnit = arguments['lengthUnit']
assert lengthUnit == 'Mb' or lengthUnit == 'gene'

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

if arguments['lengthUnit'] == 'gene':
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        nbAncGenes = len(sb.la)
        #print >> sys.stderr, "genes=%s" % nbAncGenes
        lSbsLengths.append(nbAncGenes)
elif arguments['lengthUnit'] == 'Mb':
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        lengthG1 = lengthProjectionOnGenome(1, sb, c1, genome1)
        lengthG2 = lengthProjectionOnGenome(2, sb, c2, genome2)
        averageSbLength = float(lengthG1 + lengthG2) / 2.0
        averageSbLength = averageSbLength / 1000000
        #print >> sys.stderr, "Mb=%s" % averageSbLength
        lSbsLengths.append(averageSbLength)

#X = []
#Y = []
#distribSbsLengths = collections.defaultdict(int)
#for length in lSbsLengths:
#    distribSbsLengths[length] += 1
#for (length, nbSbs) in distribSbsLengths.iteritems():
#    X.append(length)
#    Y.append(nbSbs)

#remove lengths higher than maxShownLength
lSbsLengths = [sbLength for sbLength in lSbsLengths if sbLength <= maxShownLength]
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
binwidth = 2.5
mindata=0
maxdata=80
bins = list(np.arange(mindata, maxdata + binwidth, binwidth))
n, bins, patches = ax2.hist(lSbsLengths, bins=bins, histtype='bar')
#plt.plot(X, Y, marker='o', linestyle='', color='w', label='Real value', markersize=4, markeredgewidth=1.5, markeredgecolor='k')
ax2.set_ylabel('Nb of synteny blocks')
#ax1.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], 0, "\infty"))
ax2.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], 0, maxShownLength))


# remove lengths lower than minShownLength
lSbsLengths = [sbLength for sbLength in lSbsLengths if minShownLength <= sbLength]
n, bins, patches = ax1.hist(lSbsLengths, bins=bins, histtype='bar')
t = np.asarray(bins[1:])
noisy = np.asarray(n)

# manual fit
# A=170
# K=-0.1
# C=0

# non-linear fit
#A, K, C = fit_exp_nonlinear(t, noisy)

# linear fit with the constant set to 0
C = 0
A, K = fit_exp_linear(t, noisy, C)
ysModel = model_func(t, A, K, C)
ax1.plot(bins[1:], ysModel, color='k', linewidth=1.0, label="$y = %0.2f e^{%0.2f t} + %0.2f$" % (A, K, C))
ax1.set_ylabel('Nb of synteny blocks')
ax1.set_xlabel("Lengths of synteny blocks in %s\n$%s \leq$ length $\leq %s$" % (arguments['lengthUnit'], minShownLength, maxShownLength))
ax1.legend()

#plt.tight_layout()
#plt.xlim(0, 100)
#plt.title("OSB length distribution in Human-Mouse comparison \n Confidence Interval : %s*sigma around mean" % arguments["ICfactorOfSigma"] )
#plt.title("")
fig.suptitle("Distribution of the lengths of synteny blocks\nHuman-Mouse comparison")
plt.legend()
plt.savefig(sys.stdout, format='svg')
#plt.show()
