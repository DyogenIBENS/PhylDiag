# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import math
import myTools
import operator

#############################
# Fonctions de statistiques #
#############################
class myStats:

    # Renvoie la moyenne d'une liste
    #################################
    @staticmethod
    def mean(lst):
        if len(lst) == 0:
            return None
        return float(sum(lst))/float(len(lst))

    # Renvoie l'ecart type
    #######################
    @staticmethod
    def stddev(lst, m=None, corrected=False):
        if len(lst) == 0:
            return None
        if m == None:
            m = myStats.mean(lst)
        s = sum((x-m) ** 2 for x in lst)
        return math.sqrt(float(s) / float(len(lst)-int(corrected)))

    # Renvoie la valeur comprise a x% des donnees (x=50 -> mediane)
    ################################################################
    @staticmethod
    def getValue(lst, x):
        if len(lst) == 0:
            return None
        return lst[int((x*(len(lst)-1))/100.)]

    # Renvoie la valeur telle que x% soit au dessus (x=50 -> N50)
    # Before using this function, do not forget to sort the input list
    #############################################################
    @staticmethod
    def getValueNX(lst, x):
        if len(lst) == 0:
            return None
        tmp = (sum(lst) * float(x)) / 100.
        for x in lst.__reversed__():
            tmp -= x
            if tmp <= 0:
                return x


    # Renvoie (min quart1 median quart3 max mean stddev len)
    #########################################################
    @staticmethod
    def valSummary(lst):
        l = list(lst)
        l.sort()
        m = myStats.mean(l)
        return (myStats.getValue(l, 0), myStats.getValue(l, 25),myStats.getValue(l, 50),myStats.getValue(l, 75), \
                myStats.getValueNX(l, 75),myStats.getValueNX(l, 50),myStats.getValueNX(l, 25), myStats.getValue(l, 100), m, myStats.stddev(l, m), len(l))

    @staticmethod
    def txtSummary(lst, withN50=True):
        pattern = ["%s [%s/%s/%s] "]
        if withN50:
            pattern.append("[%s/%s/%s] ")
        pattern.append("%s [")
        pattern.append("%.2f/%.2f" if len(lst) > 0 else "%s/%s")
        pattern.append("-%d]")
        if len(lst) > 0:
            data = myStats.valSummary(lst)
            return ''.join(pattern) % (data if withN50 else data[:4] + data[-4:])
        else:
            return ''.join(pattern) % ((None,)*(10 if withN50 else 7) + (0,) )


    # Renvoie la moyenne ponderee d'une liste
    ##############################
    @staticmethod
    def weightedMean(lst):
        sV = 0.
        sP = 0.
        for (val,poids) in lst:
            sV += val*poids
            sP += poids
        if sP == 0:
            return 0.
        return sV/sP

    # Renvoie la correlation entre les deux variables
    ##################################################
    @staticmethod
    def correlation(x, y):
        N = len(x)
        if N != len(y):
            N = min(N, len(y))
            x = x[:N]
            y = y[:N]

        sum_sq_x = 0.
        sum_sq_y = 0.
        sum_coproduct = 0.
        mean_x = x[0]
        mean_y = y[0]
        for i in range(1,N):
            sweep = i / (i + 1.0)
            delta_x = x[i] - mean_x
            delta_y = y[i] - mean_y
            sum_sq_x += delta_x * delta_x * sweep
            sum_sq_y += delta_y * delta_y * sweep
            sum_coproduct += delta_x * delta_y * sweep
            mean_x += delta_x / (i + 1.0)
            mean_y += delta_y / (i + 1.0)
        pop_sd_x = math.sqrt( sum_sq_x / N )
        pop_sd_y = math.sqrt( sum_sq_y / N )
        cov_x_y = sum_coproduct / N
        correlation = cov_x_y / (pop_sd_x * pop_sd_y)
        return correlation


    # Renvoie la quantite de chaque valeur dans la liste
    #####################################################
    @staticmethod
    def count(l):
        import collections
        d = collections.defaultdict(int)
        for x in l:
            d[x] += 1
        return d


######################################################
# Genere les permutations de la liste l              #
# Chaque permutation est referee par un index unique #
######################################################
class permutationGenerator:

    def __init__(self, l):

        self.l = tuple(l)
        self.n = len(l)
        fact = [1] * self.n
        for i in xrange(self.n,1,-1):
            fact[i-2] = i * fact[i-1]
        self.fact = fact
        self.a = [None] * self.n

    def getNbPerm(self, p):
        return self.fact[self.n-p-1]

    def getPermutation(self, k, p):

        assert k >= 0
        assert k < self.fact[self.n-p-1]

        for i in xrange(1,self.n):
            (self.a[i],k) = divmod(k, self.fact[i])

        b = self.l[:]
        for i in xrange(1,self.n):
            (b[i],b[self.a[i]]) = (b[self.a[i]],b[i])

        return b[-p:]


#############################
# Fonctions d'interpolation #
#############################
class myInterpolator:

    # Interpolation lineaire (1 dimension)
    ########################################
    @staticmethod
    def oneDimLinear(points, val):

        import myTools

        intervals = []

        for ((u,FU), (v,FV)) in myTools.myIterator.slidingTuple(zip(points, val)):
            A = (FU - FV) / (u - v)
            B = FU - A * u
            intervals.append((A, B))
        intervals.append(intervals[-1])

        import bisect
        def tmp(x):
            i = bisect.bisect_right(points, x)
            (A, B) = intervals[i-1]
            return (A * x + B)

        return tmp


    # Interpolation cubique (1 dimension)
    #######################################
    @staticmethod
    def oneDimCubic(points, val):

        import myTools

        intervals = []

        deriv = [0] + [(FV-FU) / (2*(v-u)) for ((u,FU), (v,FV)) in myTools.myIterator.slidingTuple(zip(points, val))] + [0]
        tangeantes = [d1+d2 for (d1,d2) in myTools.myIterator.slidingTuple(deriv)]

        # Generate a cubic spline for each interpolation interval.
        for ((u,FU,DU), (v,FV,DV)) in myTools.myIterator.slidingTuple(zip(points,val,tangeantes)):

            denom = float((u - v)**3)

            A = ((-DV - DU) * v + (DV + DU) * u + 2 * FV - 2 * FU) / denom
            B = -((-DV - 2 * DU) * v**2  + u * ((DU - DV) * v + 3 * FV - 3 * FU) + 3 * FV * v - 3 * FU * v + (2 * DV + DU) * u**2) / denom
            C = (- DU * v**3  + u * ((- 2 * DV - DU) * v**2  + 6 * FV * v - 6 * FU * v) + (DV + 2 * DU) * u**2 * v + DV * u**3) / denom
            D = -(u *(-DU * v**3  - 3 * FU * v**2) + FU * v**3 + u**2 * ((DU - DV) * v**2 + 3 * FV * v) + u**3 * (DV * v - FV)) / denom

            intervals.append((A, B, C, D))
        intervals.append(intervals[-1])

        import bisect
        def tmp(x):
            i = bisect.bisect_right(points, x)
            (A, B, C, D) = intervals[i-1]
            return ((A * x + B) * x + C) * x + D

        def c_code():
            def codeChoice(intervalList):
                n = len(intervalList)
                if n < 2:
                    return ("A=%.10e;B=%.10e;C=%.10e;D=%.10e;" % intervalList[0])
                n2 = n / 2
                return ("if (x < %.10e) {%s} else {%s}" % (points[n2], codeChoice(intervalList[:n2]), codeChoice(intervalList[n2:])))
            return ("double interpolator(double x) {double A,B,C,D;" + codeChoice(intervals) + "return ((A * x + B) * x + C) * x + D;}")

        return tmp


    # Cree un interpolateur pour chaque dimension
    ###############################################
    @staticmethod
    def getMultDim(interpolator, points, val):
        linter = [interpolator(points, [x[i] for x in val]) for i in xrange(len(val[0]))]
        return lambda x: tuple(inter(x) for inter in linter)


####################################
# Generateur de valeurs aleatoires #
####################################
class randomValue:

    import bisect
    import random

    # Renvoie un indice au hasard dans l, compte tenu des valeurs de l, qui indiquent une frequence d'apparition
    #############################################################################################################
    @staticmethod
    def bisectChooser(l):
        newl = [0] * (len(l)+1)
        x = 0
        for i in xrange(len(l)):
            x += l[i]
            newl[i+1] = newl[i] + l[i]
            assert newl[i+1] == x
        m = newl[-1]
        assert m == x
        return lambda : randomValue.bisect.bisect_left(newl, randomValue.random.random() * m) - 1


    # Distribution VonMises restreinte a [0,1] et centree sur une moyenne
    ######################################################################
    @staticmethod
    def myVonMises(mean, kappa):
        r = randomValue.random.vonmisesvariate(0, kappa) / (2*math.pi) + mean
        if r < 0:
            return mean * abs(r) / (0.5-mean)
        elif r > 1:
            return 1 - (1-mean) * abs(1-r) / (mean-0.5)
        else:
            return r

    # Tirage aleatoire selon une densite issue d'une loi geometrique
    ##################################################################
    @staticmethod
    def geometric(p):
        return int(math.ceil(math.log(1.0 - randomValue.random.random(), 1.0 - p)))


    import myTools
    @staticmethod
    @myTools.memoize
    def intParetoMean(alpha, precision, niter):
        s = 0.
        n = 0
        lastm = 0
        while True:
            for _ in xrange(niter):
                s += int(randomValue.random.paretovariate(alpha))
                n += 1
            newm = s/n
            if abs(lastm-newm) < precision:
                return newm
            lastm = newm

    @staticmethod
    def paretoAlphaFromMean(m, step=0.01, precision=0.01, niter=100000):
        alpha = int((1. + 1./(m-1.)) / step) * step
        tries = {}
        while alpha not in tries:
            res = randomValue.intParetoMean(alpha, precision, niter)
            tries[alpha] = res
            alpha += -step if res < m else step
        return min((abs(y-m),x) for (x,y) in tries.iteritems())[1]


######################
# Applatit une liste #
######################
def flatten(lst):
    res = []
    for x in lst:
        res.extend(x)
    return res


##########################################################################
# Renvoie la racine carre d'un entier, sous forme de chaine de caractere #
##########################################################################
def sqrti(n, part, nbdec):
    if (len(n) % 2) == 1:
        n = "0" + n
    backup = len(n) / 2
    iniL = len(n) + len(part)
    s = n + (part + ("00" * nbdec))[:2*nbdec]
    maxL = len(s)
    c = 0
    p = 0
    i = 0
    while ((c != 0) or (i < iniL)) and (i<maxL):
        c = 100*c + int(s[i:i+2])
        i += 2
        for x in [9,8,7,6,5,4,3,2,1,0]:
            y = (20*p+x)*x
            if y <= c:
                break
        p = 10*p+x
        c -= y
    sp = str(p)
    return (sp[:backup], sp[backup:])

################################################################
# Fast computation of combinations for probability computation #
################################################################
@myTools.memoize
def combinations(n,p):
    if n < 0:
        raise ValueError("combination(n=%s,p=%s), n cannot be < 0" % (n,p))

    if p > n:
        return 0 # there is 0 elements with p elements
    else:
        if p == 0:
            return 1 # there is only 1 combination with 0 element, the null ensembl
        elif p == n:
            return 1 # there is only 1 combination with n elements, the ensembl that conatins n elements
        else:
            try:
                # Pascal's recurrence formula, works well with the @myTools.memoize decorator
                return float(combinations(n-1,p-1) + combinations(n-1,p))
            except:
                #Pascal's recursion fails sometimes because of a "maximum recursion depth exceeded"
                num = prod([a for a in range(n-p+1,n+1)])
                den = prod([a for a in range(1,p+1)])
                try:
                    num = float(num)
                except OverflowError:
                    num = float('inf')
                try:
                    den = float(den)
                except OverflowError:
                    den = float('inf')

                if num == float('inf') and den == float('inf'):
                    return None
                elif num == float('inf') and den != float('inf'):
                    return float('inf')
                elif num != float('inf') and den == float('inf'):
                    return 0.0
                else:
                    return float(num)/den

#Warning: this function should be called from python 2.7 otherwise it may return annoying warnings in the sys.stderr as:
#Exception RuntimeError: 'maximum recursion depth exceeded in __subclasscheck__' in <type 'exceptions.RuntimeError'> ignored
@myTools.minimalPythonVersion((2,7))
def prod(factors):
    return reduce(operator.mul, factors, 1)
