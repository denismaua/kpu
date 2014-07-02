# implements SetFactor class

from math import exp, log

LOGZERO = -1000.0 #float("-inf")

class SetFactor:
    """ Class SetFactor """
    #__slots__ = ['scope','stride','dimension','zeroValue','values','id'] # not good for pypy
    def __init__(self, variableList, defaultValue=0, zeroValue=0):
        self.scope = variableList        
        self.stride = {}
        self.dimension = 1
        self.zeroValue = zeroValue
        self.defaultValue = defaultValue
        if len(self.scope) == 1:
            self.stride[self.scope[0]] = 1
            self.dimension = self.scope[0].num_states
        else:
            for variable in self.scope:
                self.stride[variable] = self.dimension
                self.dimension = self.dimension*variable.num_states
        self.num_tables = 0
        # initialize with no tables
        self.tables = []
        self.labels = []
        self.id = SetFactor.factors
        SetFactor.factors += 1

    def getStride(self,variable):
        if variable in self.stride:
            return self.stride[variable]
        else:
            return 0

    def addEmptyTable(self):
        self.tables.append( self.dimension*[self.defaultValue] )
        self.labels.append( '' )
        self.num_tables += 1
        
    def addTable(self, table):
        #assert(len(table) == self.dimension) # uncomment to assume this "unfriendly" user
        self.tables.append( table )
        self.labels.append( '' )
        self.num_tables += 1

    def removeTable(self,table):
        # TO-DO
        pass
                
    def setValue(self,table,assignment,value):
        index = 0
        for i,l in enumerate(assignment):
            index += l*self.getStride(self.scope[i])
        self.tables[table][index] = value

    def clamp(self,variable,assignment):
        # not working
        dimension = self.dimension/variable.num_states
        values = dimension*[0.0]
        j = 0
        for i in range(self.dimension):
            index = int(1.0*i/self.stride[variable]) % variable.cardinality            
            if index == assignment:
                values[j] = self.values[i]
                j += 1
            # if index != assignment:
            #    self.values[i] = self.zeroValue
        del self.stride[variable]
        self.scope.remove(variable)
        self.dimension = 1
        for variable in self.scope:
            self.stride[variable] = self.dimension
            self.dimension = self.dimension*variable.num_states
        self.values = values

    def getValue(self,table,assignment):
        index = 0
        for i,l in enumerate(assignment):
            index += l*self.getStride(self.scope[i])
        return self.tables[table][index]
    
    def printOut(self):
        print 'Dimension: %d x %d' % (self.num_tables,self.dimension)
        print 'Scope:', [var.label for var in self.scope]
        assignment = len(self.scope)*[0]
        for t,table in enumerate(self.tables):
            print "Table", t, self.labels[t]
            for i in range(self.dimension):
                for j,var in enumerate(self.scope):
                    assignment[j] = int(1.0*i/self.stride[var]) % var.num_states
                print i, assignment, table[i]

SetFactor.factors = 0

### OPERATIONS

## def SetFactorProduct(F1,F2):
##     """ Computes the product of two setfactors """
##     scope = list(set(F1.scope).union(set(F2.scope)))
##     F = SetFactor(scope)
##     for t1 in range(F1.num_tables):
##         for t2 in range(F2.num_tables):
##             F.addEmptyTable()
##             F.labels[t2+F2.num_tables*t1] = F1.labels[t1]+' '+F2.labels[t2]
##     j,k = 0,0
##     assignment = [0 for l in range(len(scope))]
##     for i in range(F.dimension):
##         for t1 in range(F1.num_tables):
##             for t2 in range(F2.num_tables):
##                 # this is very inneficient, it should be computed table by table
##                 F.tables[t2+F2.num_tables*t1][i] = F1.tables[t1][j] * F2.tables[t2][k]
##         for l in range(len(scope)):
##             assignment[l] = assignment[l] + 1
##             if assignment[l] == scope[l].num_states:
##                 assignment[l] = 0
##                 j = j - (scope[l].num_states-1)*F1.getStride(scope[l])
##                 k = k - (scope[l].num_states-1)*F2.getStride(scope[l])
##             else:
##                 j = j + F1.getStride(scope[l])
##                 k = k + F2.getStride(scope[l])
##                 break
##     return F

def SetFactorProduct(F1,F2):
    """ Computes the product of two setfactors """
    scope = list(set(F1.scope).union(set(F2.scope)))
    F = SetFactor(scope)
    j,k = 0,0
    assignment = [0 for l in range(len(scope))]
    index1 = F.dimension*[0]
    index2 = F.dimension*[0]
    for i in range(F.dimension):
        for t1 in range(F1.num_tables):
            for t2 in range(F2.num_tables):
                index1[i] = j
                index2[i] = k
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
    for t1 in range(F1.num_tables):
        for t2 in range(F2.num_tables):            
            F.addEmptyTable()
            F.labels[t2+F2.num_tables*t1] = F1.labels[t1]+F2.labels[t2]
            table = F.tables[t2+F2.num_tables*t1]
            table1 = F1.tables[t1]
            table2 = F2.tables[t2]
            for i in range(F.dimension):
                table[i] = table1[index1[i]] * table2[index2[i]]                
    return F


def SetFactorSumProduct(F1,F2,variableList):
    """ Computes the sum-product of two setfactors """
    scope = list(set(F1.scope).union(set(F2.scope)).difference(set(variableList)))
    eDim = 1
    for var in variableList:
        eDim = eDim*var.num_states        
    F = SetFactor(scope, defaultValue=0.0)
    for t1 in range(F1.num_tables):
        for t2 in range(F2.num_tables):
            F.addEmptyTable()
            F.labels[t2+F2.num_tables*t1] = F1.labels[t1]+F2.labels[t2]
    j,k = 0,0
    assignment = [0 for l in range(len(variableList)+len(scope))]
    for i in range(F.dimension):
        for _ in range(eDim):
            for t1 in range(F1.num_tables):
                for t2 in range(F2.num_tables):
                    # this is very inneficient, it should be computed table by table
                    F.tables[t2+F2.num_tables*t1][i] += F1.tables[t1][j] * F2.tables[t2][k]
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

def SetFactorSumBProduct(FactorList, variableList):
    """ Computes pairwise sum-product of list of setfactors """
    F = FactorList.pop()
    if len(FactorList) == 0:
        Fp = SetFactor([], defaultValue=1.0)
        return FactorSumProduct(F,Fp,variableList)
    
    while len(FactorList) > 1:
        Fp = FactorList.pop()
        F = FactorProduct(F,Fp)
        
    Fp = FactorList.pop()
    
    return SetFactorSumProduct(F,Fp,variableList)


def ParetoSetFactorSumBProduct(FactorList, variableList, RemoveNonExtreme=False):
    """ Computes the pareto set of pairwise sum-product of list of setfactors """
    F = FactorList.pop()
    if len(FactorList) == 0:
        Fp = SetFactor([], defaultValue=1.0)
        return ParetoSet(FactorSumProduct(F,Fp,variableList))
    
    while len(FactorList) > 1:
        Fp = FactorList.pop()
        F = ParetoSet(FactorProduct(F,Fp))
        
    Fp = FactorList.pop()
    S = ParetoSet(SetFactorSumProduct(F,Fp,variableList))
    #if RemoveNonExtreme:
    #    return ParetoSet(ExtremeSet(S))
    #if RemoveNonExtreme and S.num_tables*S.dimension < 50000:
    #    return ParetoExtremeSet(S)
    return S

def ParetoSet(F):
    """ Computes the Pareto set of the tables in a given set-factor """
    if F.num_tables == 1:
        return F

    if F.dimension == 1:
        # compute maximum value
        u = 0.0
        lbl = ''
        for i in range(F.num_tables):
            if F.tables[i][0] > u:
                u = F.tables[i][0]
                lbl = F.labels[i]
        P = SetFactor([])
        P.addTable([u])
        P.labels[0] = lbl
        return P
    
    P = SetFactor(F.scope)
    for t1 in range(F.num_tables):
        dominated = False
        for t2 in range(F.num_tables):
            if t1 != t2:
                dominated = True
                for i in range(F.dimension):
                    if F.tables[t1][i] > F.tables[t2][i]:
                        dominated = False
                        break
                if dominated:
                    break
        if not dominated:
            P.addTable(F.tables[t1])
            P.labels[-1] = F.labels[t1]
    ## numremoved= F.num_tables-P.num_tables
    ## if numremoved and P.dimension>2:
    ##     print
    ##     print numremoved, " P-dominated tables removed"
    if P.num_tables == 0:
        # all tables in F are equal
        P.addTable(F.tables[0])
        P.labels[0] = F.labels[0]
    return P

def ParetoExtremeSet(F):
    """ Compute the set of extreme non-dominated tables in a given set-factor """
    #try:
    if F.dimension == 1:
        # compute maximum value
        u = max([ F.tables[i][0] for i in range(F.num_tables) ])
        H = SetFactor([])
        H.addTable([u])
        # find a maximum strategy
        for i in range(F.num_tables):
            if F.tables[i][0] == u:
                H.labels[0] = F.labels[i]
                break
        return H
        
    if F.num_tables <= 1:
        return F
    if F.num_tables <= 10:
        return ParetoSet(F)
    
    import pulp
    H = SetFactor(F.scope)
    N = F.num_tables
    #assert(N>1)
    D = F.dimension
    for t1 in range(N):
        x = (N-1)*[None]
        for i in range(N-1):
            x[i] = pulp.LpVariable("x"+str(i), 0.0, 1.0)        
        prob = pulp.LpProblem("hull", pulp.LpMinimize)
        for i in range(D):
            c = [ F.tables[t2][i] for t2 in range(N) if t2 != t1 ]
            prob += (F.tables[t1][i] - pulp.lpDot(c,x) <= 0.0)
        prob += (pulp.lpSum(x) == 1.0)
        prob += x[0]
        status = prob.solve(pulp.GUROBI_CMD(msg=0))
        #print pulp.LpStatus[status]
        if pulp.LpStatus[status] == "Not Solved":
            # extreme point found
            H.addTable(F.tables[t1])
            H.labels[H.num_tables-1] = F.labels[t1]
    #assert( H.num_tables > 0)
    #print
    #print F.num_tables-H.num_tables, " dominated tables removed"
    return H       
    #except:
    #    print "puLP is required for linear programming. Try https://code.google.com/p/pulp-or/"
    #    exit(1)

def ExtremeSet(F):
    """ Computes the set of extreme tables in a given set-factor """
    if F.dimension == 1:
        # compute maximum value
        u = max([ F.tables[i][0] for i in range(F.num_tables) ])
        H = SetFactor([])
        H.addTable([u])
        # find a maximum strategy
        for i in range(F.num_tables):
            if F.tables[i][0] == u:
                H.labels[0] = F.labels[i]
                break
        return H
    elif F.num_tables <= F.dimension:
        return F
    else:
        import subprocess
        qhull = "qconvex" # path to qconvex
        instring = "%d\n%d\n" % (F.dimension,F.num_tables)
        for table in F.tables:
            instring += ' '.join(str(x) for x in table) + '\n'
        proc = subprocess.Popen([qhull, 'Fx', 'Pp'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        output = proc.communicate(input=instring)[0]
        if len(output) == 0:
            return F
        print
        #print F.num_tables
        #print instring
        #print
        #print output
        output = output.splitlines()
        num_vertices = int(output[0])
        assert(num_vertices > 0 and num_vertices <= F.num_tables)
        H = SetFactor(F.scope)
        for i in range(num_vertices):
            j = int(output[i+1])
            H.addTable(F.tables[j])
            H.labels[i] = F.labels[j]
        return H

def ApproxSetFactorSumBProduct(FactorList, variableList, numtables=100):
    """ Computes pairwise sum-product of list of setfactors """
    F = FactorList.pop()
    if len(FactorList) == 0:
        Fp = SetFactor([], defaultValue=1.0)
        return BucketingPruning(FactorSumProduct(F,Fp,variableList),numtables)
    
    while len(FactorList) > 1:
        Fp = FactorList.pop()
        F = BucketingPruning(FactorProduct(F,Fp),numtables)
        
    Fp = FactorList.pop()
    
    return BucketingPruning(SetFactorSumProduct(F,Fp,variableList),numtables)
    
def BucketingPruning(F, numbuckets=100):
    """ Partitions the space into backets and returns a set-potential containing at most one table per bucket """
    
    if F.num_tables <= numbuckets:
        return F

    if F.dimension == 1:
        # compute maximum value
        u = 0.0
        lbl = ""
        for i in range(F.num_tables):
            if F.tables[i][0] > u:
                u = F.tables[i][0]
                lbl = F.labels[i]
        P = SetFactor([])
        P.addTable([u])
        P.labels[0] = lbl
        return P

    # compute bucketing
    from math import pow
    scaling = int(pow(numbuckets,1.0/F.dimension))
    B = SetFactor(F.scope)
    ## buckets = {}
    ## for t,table in enumerate(F.tables):
    ##     index = ( int(table[i]*scaling) for i in range(F.dimension) )
    ##     if index not in buckets:
    ##         buckets[index] = table
    ##         B.addTable(table)
    ##         B.labels[-1] = F.labels[t]
    
    ## return B

    B = SetFactor(F.scope)
    for t1 in range(F.num_tables):
        hit = False
        for table in B.tables:
            hit = True
            for i in range(F.dimension):
                if int(F.tables[t1][i]*scaling) != int(table[i]*scaling):
                    hit = False
                    break                
            if hit:
                break
        if not hit:
            B.addTable(F.tables[t1])
            B.labels[-1] = F.labels[t1]
            if B.num_tables == numbuckets:
                break
                ## print '!!!', F.num_tables, numbuckets, B.num_tables
                ## print
                ## for table in B.tables:
                ##     print [int(table[i]*scaling) for i in range(F.dimension)], table
                ##     print [int(F.tables[t1][i]*scaling) for i in range(F.dimension)], F.tables[t1]
                ##     print
                ## print scaling
                ## exit(1)
    #assert(B.num_tables>0)
    #assert(B.num_tables <= numbuckets)
    ## numremoved= F.num_tables-P.num_tables
    ## if numremoved and P.dimension>2:
    ##     print
    ##     print numremoved, " P-dominated tables removed"
    #if P.num_tables == 0:
        # all tables in F are equal
    #    P.addTable(F.tables[0])
    return B
