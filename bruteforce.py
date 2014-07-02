# implements the multiple policy search algorithm
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
    print "[Solve ID by exhaustive search]"
    N,M=int(args[1]),int(args[2])
    # sample network
    ChanceVars, DecVars, CPT, Strategy, Utility = sampleChainID(N,M)
    run(ChanceVars, DecVars, CPT, Strategy, Utility)


def run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=True):
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
        print "trying uniform strategy..."
    stime = time.clock()
    Z, tw, wtw, mem = VariableElimination(Factors, OrderedVariables, verbose)
    etime = time.clock() - stime
    MEU = max( Z.tables[i][0] for i in range(Z.num_tables)  )
    if verbose:
        print "done. \033[91m[%gs]\033[0m" % etime
        print "Maximum memory usage: %d MB" % (mem/1000000.0)
        print "MEU:", MEU
        # save a best strategy
        for i in range(Z.num_tables):
            if Z.tables[i][0] == MEU:
                #print Z.labels[i]
                policies = Z.labels[i].split()
                for p in policies:
                    var,value = p.split('=')
                    var = int(var[1:])
                    value = int(value)
                    Strategy[var].values = Strategy[var].dimension*[0.0]
                    Strategy[var].values[value] = 1.0
                break            
    return MEU, etime
    


if __name__ == '__main__':
        if len(sys.argv) < 3:
            print "Usage:", sys.argv[0], "num_vars num_var_cardinality "
            exit(0)
 	main(sys.argv)

