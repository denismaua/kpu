# implements the multiple policy search algorithm
from SetFactor import SetFactor, SetFactorProduct, ParetoSetFactorSumBProduct, ParetoSet
from Variable import Variable
from GraphUtils import FindOrder, ComputeLowerBound, BuildDomainGraph, FindPrefix
from chainID import sampleChainID
import resource, time, sys


def VariableElimination(Factors, Ordering, verbose=True, ComputeHull=False):    
    """ Variable Elimination algorithm """
    start_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    tw = 0 # elimination width
    maxtables = 0 # max num of tables generated
    delta_mem = 0
    max_memory = 0
    cumsize = 0
    for i,var in enumerate(Ordering):
        if verbose:
        	print "-%6s\t" % var.label, 
        	sys.stdout.flush()
        B = []
        for f in Factors:
        	if var in f.scope:
        		B.append(f)
        for f in B:
        	Factors.remove(f)
    	f = ParetoSetFactorSumBProduct(B,[var], ComputeHull)
    	Factors.append(f)
        delta_mem = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) - start_mem
        max_memory = max(delta_mem, max_memory)
        maxtables = max(maxtables,f.num_tables)
        if verbose:
            dim = "%dx%d" % (f.num_tables,f.dimension)
            print "[order: %3d,\twidth: %3d,\tdim: %10s,\tsize:%10d,\tmem: %d MB,\tmaxtables: %d]" % (i,len(f.scope),dim,f.num_tables*f.dimension,delta_mem / 1000000.0,maxtables)
            sys.stdout.flush()
    f = Factors.pop()
    while len(Factors) > 0:
        fp = Factors.pop()
        f = FactorProduct(f,fp)
    return ParetoSet(f), tw, maxtables, max_memory

def main(args):
    print "[multiple policy updating]"
    N,M=int(args[1]),int(args[2])
    ChanceVars, DecVars, CPT, Strategy, Utility = sampleChainID(N,M)
    mpu,t1,numtables = run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=True, ComputeHull=True)
    print "Total elapsed time: \033[91m%gs\033[0m" % t1
    mpu,t1,numtables = run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=True, ComputeHull=False)
    print "Total elapsed time: \033[91m%gs\033[0m" % t1
    ## uncomment the following lines for chekcing correctness (for debugging on very small models only)
    ## import bruteforce
    ## e,t2 = bruteforce.run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=False)
    ## if abs(e-mpu)>0.000000001:
    ##    print e, '!=', mpu
    ##    exit(1)


def run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=True, ComputeHull=True):
    start = time.clock()
    N = len(DecVars)
    Factors = []
    # create chance setfactors (singletons)
    for n in range(len(ChanceVars)):
        f = SetFactor(CPT[n].scope)
        f.addTable(CPT[n].values)
        Factors.append(f)
    # create policy setfactors (vacuous)
    for n in range(len(DecVars)):
        f = SetFactor(Strategy[n].scope)
        for k in range(DecVars[n].cardinality):
            table = DecVars[n].cardinality*[0.0]
            table[k] = 1.0
            f.addTable(table)
            f.labels[k] = DecVars[n].label+'='+str(k)+' '
        Factors.append(f)
    # create utility setfactor
    f = SetFactor(Utility.scope)
    f.addTable(Utility.values)
    Factors.append(f)
    if verbose:
        print "building domain graph...",
        sys.stdout.flush()
    stime = time.clock()	
    dGraph = BuildDomainGraph(CPT+Strategy+[Utility], ChanceVars+DecVars)    
    etime = time.clock() - stime
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime	
        print "computing MMD lower bound on treewidth...",
        sys.stdout.flush()
    stime = time.clock()	
    mmd = ComputeLowerBound(dGraph)
    etime = time.clock() - stime
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime	
        print "applying safe reduction rules...",
        sys.stdout.flush()
    stime = time.clock()
    Prefix, low = FindPrefix(dGraph, low=mmd)
    etime = time.clock() - stime
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime	
        print "Optimal Prefix Length:", len(Prefix)
        print "computing min-fill ordering...",
        sys.stdout.flush()
    stime = time.clock()	
    OrderedVariables, tw = FindOrder(dGraph, Prefix=Prefix, treewidth=low)
    etime = time.clock() - stime
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime	
    # print "Min-Fill Elimination order:", 
    ## for var in OrderedVariables:
    ##     print var.label,
    #print
    if verbose:
        print "Treewidth:", tw
    dGraph, Prefix = None, None # dump graph
    if verbose:
        print "running variable elimination..."
    stime = time.clock()
    Z, tw, maxtables, mem = VariableElimination(Factors, OrderedVariables, verbose, ComputeHull)
    etime = time.clock() - stime
    elapsed = time.clock()-start
    assert( Z.num_tables == 1)
    MEU = Z.tables[0][0]
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime
        #print "Maximum memory usage: %d MB" % (mem/1000000.0)
        print "MEU:", MEU
        # save a best strategy
        #for i in range(Z.num_tables):
        #    if Z.tables[i][0] == MEU:
                #print Z.labels[i]
        policies = Z.labels[0].split()
        for p in policies:
            var,value = p.split('=')
            var = int(var[1:])
            value = int(value)
            Strategy[var].values = Strategy[var].dimension*[0.0]
            Strategy[var].values[value] = 1.0
#                break            
    return MEU, elapsed, maxtables
    


if __name__ == '__main__':
        if len(sys.argv) < 3:
            print "Usage:", sys.argv[0], "num_vars num_var_cardinality "
            exit(0)
 	main(sys.argv)

