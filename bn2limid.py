""" Variable Elimination Algorithm """
from Factor import Factor, LOGZERO
from Variable import Variable
from io import read_network_from_file, read_evidence_from_file
from math import exp, log

from Operations import LogFactorSumBProduct, LogFactorProduct
from Operations import FactorSumBProduct, FactorProduct


from GraphUtils import FindOrder, ComputeLowerBound, BuildDomainGraph, FindPrefix
from random import random
import resource, time

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
    print "[Covnert BN into LIMID]"
    print "loading model from file", args[1], "...",
    sys.stdout.flush()
    stime = time.clock()
    Variables, Factors = read_network_from_file(args[1], useLog=logScale)
    etime = time.clock() - stime
    print "done: %d variables, %d factors. \033[91m[%gs]\033[0m" % (len(Variables), len(Factors), etime)
    pa = {}
    ch = {}
    root = set()
    leaf = set()
    for v in Variables:
        ch[v] = set()
        pa[v] = []
    for f in Factors:
        pa[f.scope[-1]] = f.scope[:-1]
        for v in f.scope[:-1]:
            ch[v].add(f.scope[-1])    
    for v in Variables:
        if len(ch[v]) == 0:
            leaf.add(v)
        if len(pa[v]) == 0:
            root.add(v)
    print '#Leaves:', len(leaf)
    print '#Root:', len(root)
    
    print "building domain graph...",
    sys.stdout.flush()
    stime = time.clock()	
    dGraph = BuildDomainGraph(Factors, Variables)    
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
    #OrderedVariables, tw = FindOrder(dGraph)
    OrderedVariables, tw = FindOrder(dGraph, Prefix=Prefix, treewidth=low)
    etime = time.clock() - stime
    print "done. \033[91m[%gs]\033[0m" % etime
    
    print "Min-Fill Elimination order:", 
    for var in OrderedVariables:
        print var.label,
    print
    
    print "Treewidth:", tw
    dGraph, Prefix = None, None # dump graph
    prev = None
    i = 0
    for v in reversed(OrderedVariables):
        if v in leaf:
            i += 1
            #print i, ':', v.label
            w = Variable("W"+v.label[1:],2)
            f = Factor([v,w])
            util = [ random() for _ in range(v.num_states) ]
            f.values = util + [1-u for u in util]
            #f.printOut()
            Variables.append(w)
            Factors.append(f)
            o = Variable("O"+v.label[1:],2)
            Variables.append(o)
            if prev == None:
                g = Factor([w,o])
                g.values = [1.0, 0.0, 0.0, 1.0]
                Factors.append(g)
            else:
                g = Factor([w,prev,o])
                g.values = [1.0, 1.0-1.0/i, 1.0/i, 0.0, 0.0, 1.0/i, 1.0-1.0/i, 1.0]
                Factors.append(g)
            #g.printOut()
            prev = o
            #print
    V = Variable("V",1)
    U = Factor([prev,V])
    U.values = [float(i), 0.0]
    #U.printOut()
    #Factors.append(u)
    print
    # TO-DO: save to file in LIMID format and pickle dump
    out = open( args[2], "w" )
    out.write("LIMID\n" + str(len(Variables)-len(root)) + " " + str(len(root)) + " 1\n")
    v_id = {}
    i = 0
    for var in Variables:
        if var not in root:
            v_id[var.id] = i
            out.write(str(var.cardinality) + " ")
            print var.label, v_id[var.id]
            i += 1
    out.write("\n")
    for var in root:
        v_id[var.id] = i            
        out.write(str(var.cardinality) + " ")
        print var.label, v_id[var.id]
        i += 1
    out.write("\n")
    for j,var in enumerate(Variables):
        if var not in root:
            out.write(str(len(Factors[j].scope)-1) + " ")
            for pa in reversed(Factors[j].scope[:-1]):
                out.write(str(v_id[pa.id]) + " ")
            out.write("\n")
    for var in root:
        out.write("0\n")
    out.write("1 " + str(v_id[U.scope[0].id])+ "\n")

    for i,var in enumerate(Variables):
        if var not in root:
            out.write( str(len(Factors[i].values)) + "\n")
            s = Factors[i].scope
            c = len(s)*[0]
            for _ in range(Factors[i].dimension):
                out.write(" "+ str(Factors[i].getValue(c)))
                for l in reversed(range(len(c))):
                    c[l] += 1
                    if c[l] == s[l].num_states:
                        c[l] = 0
                    else:
                        break
            out.write(" \n")
    out.write( str(len(U.values)) + "\n")
    for val in U.values:
        out.write( " " + str(val))
    out.write(" \n")            
    out.close()

    
    
    ## print "running variable elimination..."	
    ## stime = time.clock()	
    ## Z, tw, wtw, mem = VariableElimination(Factors, OrderedVariables)
    ## etime = time.clock() - stime
    ## print "done. \033[91m[%gs]\033[0m" % etime
    ## Z.printOut()
    print "Elapsed time: \033[91m%gs\033[0m" % (time.clock()-start)


if __name__ == '__main__':
	import sys, time
        if len(sys.argv) < 3:
            print "Usage:", sys.argv[0], "model.UAI model.LIMID"
            exit(0)
 	main(sys.argv)
