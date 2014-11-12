# implements a k-local search on the space of policies
from SetFactor import SetFactor, SetFactorProduct, ParetoSetFactorSumBProduct, SetFactorSumBProduct, ParetoSet
from Variable import Variable
#from io import read_network_from_file, read_evidence_from_file
#from math import exp, log
from GraphUtils import FindOrder, ComputeLowerBound, BuildDomainGraph, FindPrefix
from chainID import sampleChainID
import resource, time, sys
from random import sample

def VariableElimination(Factors, Ordering, verbose=True, usePareto=True):    
    """ Variable Elimination algorithm """
    if usePareto:
        SumBProduct = ParetoSetFactorSumBProduct
    else:
        SumBProduct = SetFactorSumBProduct
    start_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    tw = 0 # elimination width
    wtw = 0 # weigthed width
    delta_mem = 0
    max_memory = 0
    max_tables = 0
    for var in Ordering:
        if verbose:
        	print "-%6s\t" % var.label, 
        	sys.stdout.flush()
        B = []
        for f in Factors:
        	if var in f.scope:
        		B.append(f)
        for f in B:
        	Factors.remove(f)
            
    	f = SumBProduct(B,[var])
        max_tables = max(max_tables,f.num_tables)
    	Factors.append(f)
        delta_mem = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) - start_mem
        max_memory = max(delta_mem, max_memory)
        if verbose:
            dim = "%dx%d" % (f.num_tables,f.dimension)
            print "[width: %3d,\tdim: %10s,\tsize:%10d,\tmem: %d MB]" % (len(f.scope),dim,f.num_tables*f.dimension,delta_mem / 1000000.0)
            tw = max(tw, len(f.scope))
            wtw = max(wtw, f.dimension)
            sys.stdout.flush()
    f = Factors.pop()
    while len(Factors) > 0:
        fp = Factors.pop()
        f = SetFactorProduct(f,fp)
        max_tables = max(max_tables,f.num_tables)
    return ParetoSet(f), tw, wtw, max_memory, max_tables

def main(args):
    print "[Solve ID by local search]"
    N,M,K=int(args[1]),int(args[2]),int(args[3])
    from Factor import Factor
    # sample network
    ChanceVars, DecVars, CPT, Strategy, Utility = sampleChainID(N,M)
    #print "loading model from file...",
    #stime = time.clock()
    #Variables, Factors = read_network_from_file(args[1], useLog=logScale)
    #etime = time.clock() - stime
    #print "done: %d variables, %d factors. \033[91m[%gs]\033[0m" % (len(Variables), len(Factors), etime)
    run(ChanceVars, DecVars, CPT, Strategy, Utility, K, False, 100, True)

    sol = len(DecVars)*[ None ]
    for n in range(len(DecVars)):
        sol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
    kPUvalue, kPUit, kPUtime, KPUsize = run(ChanceVars, DecVars, CPT, sol, Utility, K, False, 100, True)
    sol2 = len(DecVars)*[ None ]
    for n in range(len(DecVars)):
        sol2[n] = Factor([DecVars[n]], defaultValue=1.0/M)
    kPUvalue2, kPUit2, kPUtime2, KPUsize2 = run(ChanceVars, DecVars, CPT, sol2, Utility, K, False, 100, False)
    print
    print "METHOD    \tTIME(s)   \tMAX SIZE  \tVALUE"
    print "Pruned    \t%10s\t%10s\t%g" % (str(kPUtime), str(KPUsize), kPUvalue)
    print "Exhaustive\t%10s\t%10s\t%g" % (str(kPUtime2), str(KPUsize2),  kPUvalue2)

    
    
def run(ChanceVars, DecVars, CPT, Strategy, Utility, NSize, verbose=True, maxiterations=50, usePareto=True):
    start = time.clock()
    if verbose:
        sys.stdout.flush()
    N = len(Strategy)
    ## if verbose:
    ##     print "building domain graph...",
    ##     sys.stdout.flush()
    ## stime = time.clock()	
    ## dGraph = BuildDomainGraph(CPT+Strategy+[Utility], ChanceVars+DecVars)    
    ## etime = time.clock() - stime
    ## if verbose:
    ##     print "done. \033[91m[%gs]\033[0m" % etime	
    ##     print "computing MMD lower bound on treewidth...",
    ##     sys.stdout.flush()
    ## stime = time.clock()	
    ## mmd = ComputeLowerBound(dGraph)
    ## etime = time.clock() - stime
    ## if verbose:
    ##     print "done. \033[91m[%gs]\033[0m" % etime	
    ##     print "applying safe reduction rules...",
    ##     sys.stdout.flush()
    ## stime = time.clock()
    ## Prefix, low = FindPrefix(dGraph, low=mmd)
    ## etime = time.clock() - stime
    ## if verbose:
    ##     print "done. \033[91m[%gs]\033[0m" % etime	
    ##     print "Optimal Prefix Length:", len(Prefix)
    ##     print "computing min-fill ordering...",
    ##     sys.stdout.flush()
    ## stime = time.clock()	
    ## OrderedVariables, tw = FindOrder(dGraph, Prefix=Prefix, treewidth=low)
    ## etime = time.clock() - stime
    ## if verbose:
    ##     print "done. \033[91m[%gs]\033[0m" % etime	
    ## # print "Min-Fill Elimination order:", 
    ## ## for var in OrderedVariables:
    ## ##     print var.label,
    ## #print
    ## if verbose:
    ##     print "Treewidth:", tw
    ## dGraph, Prefix = None, None # dump graph
    OrderedVariables = DecVars + ChanceVars[::-1] # reverse topological order

    # create chance setfactors
    CFactors = []
    for n in range(len(ChanceVars)):
        f = SetFactor(CPT[n].scope)
        f.addTable(CPT[n].values)
        CFactors.append(f)
    # create policy setfactors (uniform)
    DFactors = []
    for n in range(len(DecVars)):
        f = SetFactor(Strategy[n].scope,defaultValue=1.0/DecVars[n].cardinality)
        f.addTable(Strategy[n].values) # uncomment to allow intial point to be passed as argument
        #f.addEmptyTable() # uncomment to for uniform as initial strategy
        DFactors.append(f)
    # create utility setfactor
    UFactor = SetFactor(Utility.scope)
    UFactor.addTable(Utility.values)

    # baseline value (clearly non-optimal)
    if verbose:
        print "trying initial strategy..."
    stime = time.clock()	
    Z, tw, wtw, mem, max_tables = VariableElimination(CFactors+DFactors+[UFactor], OrderedVariables, False)
    etime = time.clock() - stime
    MEU = Z.tables[0][0]
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime
        print "Maximum memory usage: %d MB" % (mem/1000000.0)
        print "Incumbent:", MEU
    
    # maximum number of iterations
    iteration = 0
    previousMEU = 0.0
    while iteration < maxiterations: # and previousMEU != MEU:
        previousMEU = MEU
        if verbose:
            print "#%d" % (iteration + 1)

        # sample k decision variables
        d = iteration % N
        others = [i for i in xrange(N) if i != d]
        ksearch = sample(others,NSize-1)+[d] # ensure all variables are selected at least once (with enough iterations)
        # find best action for D[ksearch]
        if verbose:
            print "optimizing for", ' '.join(DecVars[d].label for d in ksearch), "..."
        Policies = []
        for d in ksearch:
            P = SetFactor( [ DecVars[d] ], defaultValue=0.0 )
            for j in range(DecVars[d].cardinality):
                P.addEmptyTable()
                P.tables[j][j] = 1.0
                P.labels[j] = ' '+DecVars[d].label+'='+str(j)
            Policies.append(P)
        StaticPolicies = [ DFactors[n] for n in range(N) if n not in ksearch ]
        # check expected value of new strategy
        if verbose:
            print "running variable elimination..."
        stime = time.clock()	
        Z, tw, wtw, mem, max_tables = VariableElimination(CFactors+StaticPolicies+Policies+[UFactor], OrderedVariables, False, usePareto)
        etime = time.clock() - stime
        assert(Z.num_tables == 1)
        E = Z.tables[0][0] 
        if verbose:
            print "done. \033[91m[%gs]\033[0m" % etime
            print "Max Tables: %d" % max_tables
            print "Memory usage: %d MB" % (mem/1000000.0)
            ## print "Effective weigthed treewidth:", wtw
            print "Improvement:", (E-MEU)
        if E > MEU:
            # save best policies
            policies = Z.labels[0].split()
            #print Z.labels[0]
            for p in policies:
                var,value = p.split('=')
                var = int(var[1:])
                assert(var in ksearch)
                value = int(value)
                for i in range(Strategy[var].dimension):
                    Strategy[var].values[i] = 0.0
                Strategy[var].values[value] = 1.0
            
            MEU = E
                            
        iteration += 1
        if verbose:
            print "Incumbent:", MEU
    
    totaltime = time.clock()-start
    if verbose:
        print "Total elapsed time: \033[91m%gs\033[0m" % totaltime
        print "Iterations:", iteration        
        print "Best solution:", MEU
    return MEU, iteration, totaltime, max_tables
    
if __name__ == '__main__':
        if len(sys.argv) < 4:
            print "Usage:", sys.argv[0], "num_vars num_var_cardinality neighborhood-size"
            exit(0)
 	main(sys.argv)
