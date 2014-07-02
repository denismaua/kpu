""" Graph Utilities """
from Factor import Factor, LOGZERO

def BuildDomainGraph(Factors, Variables):
    """ Returns the domaing graph of the graphical model """
    dGraph = {}     #  domain graph (represented as an adjancency list)
    for var in Variables:
        dGraph[var] = []

    for f in Factors:
        for i in range(len(f.scope)-1):
            #if f.scope[i] in dGraph:
                for j in range(i+1,len(f.scope)):
                    #if f.scope[j] in dGraph:
                        if f.scope[j] not in dGraph[f.scope[i]]:
                            dGraph[f.scope[i]].append(f.scope[j])
                        if f.scope[i] not in dGraph[f.scope[j]]:
                            dGraph[f.scope[j]].append(f.scope[i])
    return dGraph

def SafeRules(dGraph, Painted, low):
    """ Applies safe reduction rules """
    for var in Painted:
        if not Painted[var]:
            degree = len(set(ne for ne in dGraph[var] if not Painted[ne])) # node degree in the current graph
            if degree < 2:
                # Islet or Twig 
                return var, low
            if degree < 4 and low >= degree:
                # Series or Triangle
                return var, low
            
            # Simplicial
            nolink = False
            for n1 in range(len(dGraph[var])-1):
                if not Painted[dGraph[var][n1]]:
                    for n2 in range(n1+1,len(dGraph[var])):
                        if not Painted[dGraph[var][n2]]:
                            if dGraph[var][n2] not in dGraph[dGraph[var][n1]]:
                                # n1 and n2 are not connected
                                nolink = True
                                break
                    if nolink:
                        break
            else:
                # Simplicial
                return var, max(degree, low)
                                
            # Almost Simplicial
            if low >= degree:
                for ne in range(len(dGraph[var])):
                    if not Painted[dGraph[var][ne]]:
                        # verify whether remaining neighbors form a clique
                        nolink = False
                        for n1 in range(len(dGraph[var])-1):
                            if n1 != ne and not Painted[dGraph[var][n1]]:
                                for n2 in range(n1+1,len(dGraph[var])):
                                    if n2 != ne and not Painted[dGraph[var][n2]]:
                                        if dGraph[var][n2] not in dGraph[dGraph[var][n1]]:
                                            # n1 and n2 are not connected
                                            nolink = True
                                            break
                                if nolink:
                                    break
                        else:
                            # Almost simplicial
                            return var, low
    return None, low


def MinFill(dGraph, Painted):
    """ Finds best node to eliminate according to min-fill heuristic """
    bestScore = len(dGraph)
    best = None
    for var in Painted:
        if not Painted[var]:
            score = 0
            for n1 in range(len(dGraph[var])-1):
                if not Painted[dGraph[var][n1]]:
                    for n2 in range(n1+1,len(dGraph[var])):
                        if not Painted[dGraph[var][n2]]:
                            if dGraph[var][n2] not in dGraph[dGraph[var][n1]]:
                                # n2 is not a neighbor of n1
                                score += 1
            if score < bestScore:
                bestScore = score
                best = var
    return best, bestScore

def MinDegree(dGraph, Painted):
    """ Finds best node to eliminate according to min-degree heuristic """
    bestScore = len(dGraph)
    best = None
    # for var in Unnordered:
    for var in Painted:
        if not Painted[var]:
            score = len(set(ne for ne in dGraph[var] if not Painted[ne])) # node degree in the current graph
            if score < bestScore:
                bestScore = score
                best = var
    return best, bestScore
            

def FindPrefix(dGraph, low=1):
    """ Applies safe reduction rules to generate an optimal prefix """
    Ordered = []
    Painted = dict.fromkeys(dGraph.keys(), False)
    treewidth = low
    while len(Ordered) < len(dGraph):
        # look for simplicial nodes
        best, low = SafeRules(dGraph, Painted, low)
        if best is None:
            # no safe rules to apply
            return Ordered, treewidth
        treewidth = max( treewidth, low )
        # add variable to ordered list
        Ordered.append( best )
        Painted[best] = True        
        degree = len(set(var for var in dGraph[best] if not Painted[var])) # node degree in the current graph
        treewidth = max( treewidth, degree )
        # connect unnorderd neighbors
        for n1 in range(len(dGraph[best])-1):
            if not Painted[dGraph[best][n1]]:
                for n2 in range(n1+1,len(dGraph[best])):
                    if not Painted[dGraph[best][n2]]:
                        if dGraph[best][n2] not in dGraph[dGraph[best][n1]]:
                            # n2 is not a neighbor of n1
                            dGraph[dGraph[best][n1]].append(dGraph[best][n2])
                            dGraph[dGraph[best][n2]].append(dGraph[best][n1])
    
    return Ordered, treewidth


def FindOrder(dGraph, Heuristic = MinFill, Prefix = [], treewidth=1):
    """ Computes an elimination order """
    Ordered = len(dGraph)*[None]
    Painted = dict.fromkeys(dGraph.keys(), False)
    order = 0
    # mark precomputed prefix
    for var in Prefix:
        Ordered[order] = var
        Painted[var] = True
        order += 1

    while order < len(dGraph):
        # select variable to eliminate
        best, score = Heuristic(dGraph, Painted) # find best scoring variable   
        # add variable to ordered list
        Ordered[order] = best
        order += 1
        Painted[best] = True        
        degree = len(set(var for var in dGraph[best] if not Painted[var])) # node degree in the current graph
        treewidth = max( treewidth, degree )
        # connect unnorderd neighbors
        for n1 in range(len(dGraph[best])-1):
            if not Painted[dGraph[best][n1]]:
                for n2 in range(n1+1,len(dGraph[best])):
                    if not Painted[dGraph[best][n2]]:
                        if dGraph[best][n2] not in dGraph[dGraph[best][n1]]:
                            # n2 is not a neighbor of n1
                            dGraph[dGraph[best][n1]].append(dGraph[best][n2])
                            dGraph[dGraph[best][n2]].append(dGraph[best][n1])
    
    return Ordered, treewidth    

def ComputeLowerBound(dGraph, Heuristic = MinDegree):
    """ Computes an elimination order """
    Painted = dict.fromkeys(dGraph.keys(), False)
    treewidth = 0
    for node in dGraph:
    # while len(Unnordered) > 0:
        # select variable to eliminate
        best, score = Heuristic(dGraph, Painted) # find best scoring variable   
        # add variable to ordered list
        Painted[best] = True        
        degree = len(set(var for var in dGraph[best] if not Painted[var])) # node degree in the current graph
        treewidth = max( treewidth, degree )
    
    return treewidth    


def ComputeWidth(dGraph, Order):
    """ Computes the width of an elimination order """
    Painted = dict.fromkeys(dGraph.keys(), False)
    treewidth = 0
    for var in Order:
        Painted[var] = True        
        degree = len(set(x for x in dGraph[var] if not Painted[x])) # node degree in the current graph
        treewidth = max( treewidth, degree )
        # connect unnorderd neighbors
        for n1 in range(len(dGraph[var])-1):
            if not Painted[dGraph[var][n1]]:
                for n2 in range(n1+1,len(dGraph[var])):
                    if not Painted[dGraph[var][n2]]:
                        if dGraph[var][n2] not in dGraph[dGraph[var][n1]]:
                            # n2 is not a neighbor of n1
                            dGraph[dGraph[var][n1]].append(dGraph[var][n2])
                            dGraph[dGraph[var][n2]].append(dGraph[var][n1])
    
    return treewidth    


def Topological(DAG):
    """ Computes a topological ordering of nodes in a DAG represented as a parenthood list """
    # Builds chidrenOf relation
    CH = {} # adjancecy list
    for v in DAG:
        # initialize list
        CH[v] = []
    for v in DAG:
        # find childrenOf relations
        for p in DAG[v]:
            CH[p].append(v)
    # implements Kahn's algorithm
    D = DAG.copy() # copy it so we can modify it
    # 1. Find root nodes
    L = []
    S = []
    for v in D:
        if not len(D[v]):
            S.append(v)

    while len(S):
        s = S.pop(0)
        L.append(s)
        for v in CH[s]:
            D[v].remove(s)
            if not len(D[v]):
                S.append(v)
        del CH[s]

    return L
            
