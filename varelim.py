""" Variable Elimination Algorithm """
from Factor import Factor, LOGZERO
from Variable import Variable
from io import read_network_from_file, read_evidence_from_file
from math import exp, log
from Operations import LogFactorSumBProduct, LogFactorProduct
from Operations import FactorSumBProduct, FactorProduct
from GraphUtils import FindOrder, ComputeLowerBound, BuildDomainGraph, FindPrefix
import resource, time

def LogVariableElimination(Factors, Ordering, verbose=True):    
    """ Variable Elimination algorithm """
    start_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    tw = 0 # elimination width
    wtw = 0 # weigthed width
    delta_mem = 0
    max_memory = 0
    for var in Ordering:
        if verbose:
        	print "-%s" % var.label, 
        	sys.stdout.flush()
        B = []
        for f in Factors:
        	if var in f.scope:
        		B.append(f)
        for f in B:
        	Factors.remove(f)
    	f = LogFactorSumBProduct(B,[var])
    	Factors.append(f)
        delta_mem = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) - start_mem
        max_memory = max(delta_mem, max_memory)
        if verbose:
            print "[tw: %d,\tdim: %d,\tmem: %d MB]" % (len(f.scope),f.dimension,delta_mem / 1000000.0)
            tw = max(tw, len(f.scope))
            wtw = max(wtw, f.dimension)
            sys.stdout.flush()
    f = Factors.pop()
    while len(Factors) > 0:
        fp = Factors.pop()
        f = LogFactorProduct(f,fp)
    if verbose:
        print
    return f, tw, wtw, max_memory

def VariableElimination(Factors, Ordering, verbose=True):    
    """ Variable Elimination algorithm """
    start_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if verbose:
        import sys
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
        # if len(B) < 1:
        #     print
        #     for f in Factors:
        #         print f.id, [v.label for v in f.scope]
        #     exit(1)
                
        for f in B:
        	Factors.remove(f)
    	f = FactorSumBProduct(B,[var])
    	Factors.append(f)
        delta_mem = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) - start_mem
        max_memory = max(delta_mem, max_memory)
        if verbose:
            print "[tw: %3d,\tdim: %10d,\tmem: %d MB]" % (len(f.scope),f.dimension,delta_mem / 1000000.0)
            tw = max(tw, len(f.scope))
            wtw = max(wtw, f.dimension)
            sys.stdout.flush()
    f = Factors.pop()
    while len(Factors) > 0:
        fp = Factors.pop()
        f = FactorProduct(f,fp)
    #if verbose:
        #print
    return f, tw, wtw, max_memory


def main(args):
    logScale = False # Use log-scale
    #print "[Compute all marginals by Variable Elimination]"
    start = time.clock()
    print "[Compute log-probability of evidence]"
    print "loading model from file...",
    sys.stdout.flush()
    stime = time.clock()
    Variables, Factors = read_network_from_file(args[1], useLog=logScale)
    etime = time.clock() - stime
    print "done: %d variables, %d factors. \033[91m[%gs]\033[0m" % (len(Variables), len(Factors), etime)
    if len(args) > 2:
            stime = time.clock()
            print "loading evidence from file...",
            sys.stdout.flush()
            Evidence = read_evidence_from_file(args[2])
            etime = time.clock() - stime		
            print "done: %d observations. \033[91m[%gs]\033[0m" % (len(Evidence), etime)
            print "conditioning factors...",
            sys.stdout.flush()
            stime = time.clock()
            for e in Evidence:
                    #print Variables[e].label,":=",Evidence[e]
                    for f in Factors:
                            if Variables[e] in f.scope:
                                    # clamp variable
                                    f.clamp(Variables[e], Evidence[e])
            etime = time.clock() - stime		
            print "done. \033[91m[%gs]\033[0m" % etime
    else:
        Evidence = []
    NonEvidenceVariables = [ Variables[i] for i in range(len(Variables)) if i not in Evidence ]
    print "building domain graph...",
    sys.stdout.flush()
    stime = time.clock()	
    dGraph = BuildDomainGraph(Factors, NonEvidenceVariables)    
    etime = time.clock() - stime
    print "done. \033[91m[%gs]\033[0m" % etime	
    print "computing MMD lower bound on treewidth...",
    sys.stdout.flush()
    stime = time.clock()	
    mmd = ComputeLowerBound(dGraph)
    etime = time.clock() - stime
    print "done. \033[91m[%gs]\033[0m" % etime	
    print "applying safe reduction rules...",
    sys.stdout.flush()
    stime = time.clock()
    Prefix, low = FindPrefix(dGraph, low=mmd)
    etime = time.clock() - stime
    print "done. \033[91m[%gs]\033[0m" % etime	
    print "Optimal Prefix Length:", len(Prefix)
    print "computing min-fill ordering...",
    sys.stdout.flush()
    stime = time.clock()	
    OrderedVariables, tw = FindOrder(dGraph, Prefix=Prefix, treewidth=low)
    etime = time.clock() - stime
    print "done. \033[91m[%gs]\033[0m" % etime	
    # print "Min-Fill Elimination order:", 
    # for var in OrderedVariables:
    #     print var.label,
    # print
    print "Treewidth:", tw
    dGraph, Prefix = None, None # dump graph
    print "running variable elimination..."	
    stime = time.clock()	
    if logScale:
        Z, tw, wtw, mem = LogVariableElimination(Factors, OrderedVariables)
    else:
        Z, tw, wtw, mem = VariableElimination(Factors, OrderedVariables)
    etime = time.clock() - stime
    print "done. \033[91m[%gs]\033[0m" % etime	
    print "Maximum Memory Usage: %d MB" % (mem/1000000.0)
    print "Effective treewidth:", tw
    print "Effective weigthed treewidth:", wtw
    if logScale:
        print "\033[94mlog_10 Z: %g\033[0m" % (Z.values[0]/log(10)) # partition function is defined as log_10 Pr(evidence) 
    else:
        print "\033[94mlog_10 Z: %g\033[0m" % (log(Z.values[0])/log(10)) # partition function is defined as log_10 Pr(evidence) 
    #print "\033[94mP(evidence): %g" % exp(Z.values[0])
    # query = Variables[36]
    # M = VariableElimination(Factors, [ var for var in OrderedVariables if var != query ])
    # for i in range(M.dimension):
    # 	M.values[i] = exp(M.values[i])
    # Z = sum(M.values)	
    # print "\033[94mP(%s|evidence):" % query.label, # posterior probability
    # for i in range(M.dimension):
    # 	print M.values[i]/Z,
    # print "\033[0m"
    # #M.printOut()
    print "Elapsed time: \033[91m%gs\033[0m" % (time.clock()-start)


if __name__ == '__main__':
	import sys, time
        if len(sys.argv) < 2:
            print "Usage:", sys.argv[0], "model.UAI [evidence.UAI]"
            exit(0)
 	main(sys.argv)
