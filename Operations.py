# implements basic operations with tabular factors

from Factor import Factor, LOGZERO
from Variable import Variable
from math import exp, log

def FactorProduct(F1,F2):
    """ Computes the product of two factors """
    scope = list(set(F1.scope).union(set(F2.scope)))
    F = Factor(scope)
    j,k = 0,0
    assignment = [0 for l in range(len(scope))]
    for i in range(F.dimension):
        F.values[i] = F1.values[j] * F2.values[k]
        for l in range(len(scope)):
            assignment[l] = assignment[l] + 1
            if assignment[l] == scope[l].num_states:
                assignment[l] = 0
                j = j - (scope[l].num_states-1)*F1.getStride(scope[l])
                k = k - (scope[l].num_states-1)*F2.getStride(scope[l])
            else:
                j = j + F1.getStride(scope[l])
                k = k + F2.getStride(scope[l])
                break
    return F

def LogFactorProduct(F1,F2):
    """ Computes the product of two factors specified in log scale """
    scope = list(set(F1.scope).union(set(F2.scope)))
    F = Factor(scope)
    j,k = 0,0
    assignment = [0 for l in range(len(scope))]
    for i in range(F.dimension):
        F.values[i] = F1.values[j] + F2.values[k]
        for l in range(len(scope)):
            assignment[l] = assignment[l] + 1
            if assignment[l] == scope[l].num_states:
                assignment[l] = 0
                j = j - (scope[l].num_states-1)*F1.getStride(scope[l])
                k = k - (scope[l].num_states-1)*F2.getStride(scope[l])
            else:
                j = j + F1.getStride(scope[l])
                k = k + F2.getStride(scope[l])
                break
    return F

def FactorSumProduct(F1,F2,variableList):
    """ Computes the sum-product of two factors w.r.t. the variables in variableList"""
    scope = list(set(F1.scope).union(set(F2.scope)).difference(set(variableList)))
    eDim = 1
    for var in variableList:
        eDim = eDim*var.num_states        
    F = Factor(scope, defaultValue=0.0)
    j,k = 0,0
    assignment = [0 for l in range(len(variableList)+len(scope))]
    for i in range(F.dimension):
        for _ in range(eDim):
            F.values[i] += F1.values[j] * F2.values[k]
            for l in range(len(variableList)+len(scope)):
                assignment[l] = assignment[l] + 1
                if l < len(variableList):
                    var = variableList[l]
                else:
                    var = scope[l-len(variableList)]
                if assignment[l] == var.num_states:
                    assignment[l] = 0
                    j = j - (var.num_states-1)*F1.getStride(var)
                    k = k - (var.num_states-1)*F2.getStride(var)
                else:
                    j = j + F1.getStride(var)
                    k = k + F2.getStride(var)
                    break
    return F

def LogFactorSumProduct(F1,F2,variableList):
    """ Computes the sum-product of two factors specified in log-scale """
    scope = list(set(F1.scope).union(set(F2.scope)).difference(set(variableList)))
    eDim = 1
    for var in variableList:
        eDim = eDim*var.num_states        
    F = Factor(scope, defaultValue=0.0)
    j,k = 0,0
    c1 = max(F1.values)
    c2 = max(F2.values)
    c = max(c1,c2)
    assignment = [0 for l in range(len(variableList)+len(scope))]
    for i in range(F.dimension):
        for _ in range(eDim):
            F.values[i] += exp(F1.values[j]+F2.values[k]-c)
            for l in range(len(variableList)+len(scope)):
                assignment[l] = assignment[l] + 1
                if l < len(variableList):
                    var = variableList[l]
                else:
                    var = scope[l-len(variableList)]
                if assignment[l] == var.num_states:
                    assignment[l] = 0
                    j = j - (var.num_states-1)*F1.getStride(var)
                    k = k - (var.num_states-1)*F2.getStride(var)
                else:
                    j = j + F1.getStride(var)
                    k = k + F2.getStride(var)
                    break
    for k in range(F.dimension):
        try:
        	F.values[k] = log(F.values[k]) + c
        except ValueError:
        	F.values[k] = LOGZERO
    return F

def FactorMaxProduct(F1,F2,variableList):
    """ Computes the max-product of two factors specified """
    scope = list(set(F1.scope).union(set(F2.scope)).difference(set(variableList)))
    eDim = 1
    for var in variableList:
        eDim = eDim*var.num_states        
    F = Factor(scope, defaultValue=0.0)
    j,k = 0,0
    assignment = [0 for l in range(len(variableList)+len(scope))]
    for i in range(F.dimension):
        for _ in range(eDim):
            F.values[i] = max(F.values[i], F1.values[j]*F2.values[k])
            for l in range(len(variableList)+len(scope)):
                assignment[l] = assignment[l] + 1
                if l < len(variableList):
                    var = variableList[l]
                else:
                    var = scope[l-len(variableList)]
                if assignment[l] == var.num_states:
                    assignment[l] = 0
                    j = j - (var.num_states-1)*F1.getStride(var)
                    k = k - (var.num_states-1)*F2.getStride(var)
                else:
                    j = j + F1.getStride(var)
                    k = k + F2.getStride(var)
                    break
    return F

def LogFactorMaxProduct(F1,F2,variableList):
    """ Computes the max-product of two factors specified in log-scale """
    scope = list(set(F1.scope).union(set(F2.scope)).difference(set(variableList)))
    eDim = 1
    for var in variableList:
        eDim = eDim*var.num_states        
    F = Factor(scope)
    j,k = 0,0
    assignment = [0 for l in range(len(variableList)+len(scope))]
    for i in range(F.dimension):
        for _ in range(eDim):
            F.values[i] = max(F.values[i], F1.values[j]+F2.values[k])
            for l in range(len(variableList)+len(scope)):
                assignment[l] = assignment[l] + 1
                if l < len(variableList):
                    var = variableList[l]
                else:
                    var = scope[l-len(variableList)]
                if assignment[l] == var.num_states:
                    assignment[l] = 0
                    j = j - (var.num_states-1)*F1.getStride(var)
                    k = k - (var.num_states-1)*F2.getStride(var)
                else:
                    j = j + F1.getStride(var)
                    k = k + F2.getStride(var)
                    break
    return F

def FactorSumBProduct(FactorList, variableList):
    """ Computes pairwise sum-product of list of factors """
    F = FactorList.pop()
    if len(FactorList) == 0:
        Fp = Factor([], defaultValue=1.0)
        return FactorSumProduct(F,Fp,variableList)
    
    while len(FactorList) > 1:
        Fp = FactorList.pop()
        F = FactorProduct(F,Fp)

    Fp = FactorList.pop()
    
    return FactorSumProduct(F,Fp,variableList)

def LogFactorSumBProduct(FactorList, variableList):
    """ Computes pairwise sum-product of list of factors specified in log-scale """
    F = FactorList.pop()
    if len(FactorList) == 0:
        Fp = Factor([], defaultValue=0.0)
        return LogFactorSumProduct(F,Fp,variableList)
    
    while len(FactorList) > 1:
        Fp = FactorList.pop()
        F = LogFactorProduct(F,Fp)

    Fp = FactorList.pop()
    
    return LogFactorSumProduct(F,Fp,variableList)

def FactorMaxBProduct(FactorList, variableList):
    """ Computes pairwise max-product of list of factors """
    F = FactorList.pop()
    if len(FactorList) == 0:
        Fp = Factor([], defaultValue=1.0)
        return FactorMaxProduct(F,Fp,variableList)
    
    while len(FactorList) > 1:
        Fp = FactorList.pop()
        F = FactorProduct(F,Fp)

    Fp = FactorList.pop()
    
    return FactorMaxProduct(F,Fp,variableList)
