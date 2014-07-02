# implements a local search on the space of policies
from SetFactor import SetFactor, SetFactorProduct, SetFactorSumBProduct
from Variable import Variable
#from io import read_network_from_file, read_evidence_from_file
#from math import exp, log
from GraphUtils import FindOrder, ComputeLowerBound, BuildDomainGraph, FindPrefix
from chainID import sampleChainID
import resource, time, sys

def VariableElimination(Factors, Ordering, verbose=True):    
    """ Variable Elimination algorithm """
    start_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    tw = 0 # elimination width
    wtw = 0 # weigthed width
    delta_mem = 0
    max_memory = 0
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
    	f = SetFactorSumBProduct(B,[var])
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
    return f, tw, wtw, max_memory

def main(args):
    print "[Solve ID by local search]"
    N,M=int(args[1]),int(args[2])
    # sample network
    ChanceVars, DecVars, CPT, Strategy, Utility = sampleChainID(N,M)
    #print "loading model from file...",
    #stime = time.clock()
    #Variables, Factors = read_network_from_file(args[1], useLog=logScale)
    #etime = time.clock() - stime
    #print "done: %d variables, %d factors. \033[91m[%gs]\033[0m" % (len(Variables), len(Factors), etime)
    MEU, iteration, time = run(ChanceVars, DecVars, CPT, Strategy, Utility)
    
def run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=True):
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
    Z, tw, wtw, mem = VariableElimination(CFactors+DFactors+[UFactor], OrderedVariables, False)
    etime = time.clock() - stime
    MEU = Z.tables[0][0]
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime
        print "Memory usage: %d MB" % (mem/1000000.0)
        print "Incumbent:", MEU
    
    # maximum number of iterations
    maxiterations = 100
    iteration = 0
    previousMEU = 0.0
    while iteration < maxiterations and previousMEU != MEU:
        previousMEU = MEU
        if verbose:
            print "#%d" % (iteration + 1)
        for i in range(N):
            d = N-i-1
            # find best action for D[d]
            if verbose:
                print "optimizing for", DecVars[d].label,"...",
                sys.stdout.flush()
            Policies = SetFactor( [ DecVars[d] ] )
            for j in range(DecVars[d].cardinality):
                Policies.addEmptyTable()
                Policies.tables[j][j] = 1.0
                Policies.labels[j] = ' '+DecVars[d].label+'='+str(j)
            # check expected value of new strategy
            #if verbose:
            #    print "running variable elimination..."
            stime = time.clock()	
            Z, tw, wtw, mem = VariableElimination(CFactors+DFactors[:d]+[Policies]+DFactors[d+1:]+[UFactor], OrderedVariables, False)
            etime = time.clock() - stime
            E = max( Z.tables[i][0] for i in range(Z.num_tables)  )
            if verbose:
                print "done. \033[91m[%gs]\033[0m" % etime
                print "Memory usage: %d MB" % (mem/1000000.0)
                ## print "Effective treewidth:", tw
                ## print "Effective weigthed treewidth:", wtw
                print "Improvement:", E-MEU

            if E > MEU:
                for i in range(Z.num_tables):
                    if Z.tables[i][0] == E:
                        break
                bestAction = int(Z.labels[i].split('=')[1])
                # save a best policy for D[d]
                for k in range(DecVars[d].cardinality):
                    if k == bestAction:
                        #DFactors[d].tables[0][k] = 1.0
                        Strategy[d].values[k] = 1.0
                    else:
                        #DFactors[d].tables[0][k] = 0.0
                        Strategy[d].values[k] = 0.0
                MEU = E

            #Z.printOut()
                            
        iteration += 1
        if verbose:
            print "Incumbent:", MEU        
    totaltime = time.clock()-start
    if verbose:
        print "Total elapsed time: \033[91m%gs\033[0m" % totaltime
        print "Iterations:", iteration        
        print "Best solution:", MEU
    return MEU, iteration, totaltime
    
if __name__ == '__main__':
        if len(sys.argv) < 3:
            print "Usage:", sys.argv[0], "num_vars num_var_cardinality "
            exit(0)
 	main(sys.argv)
