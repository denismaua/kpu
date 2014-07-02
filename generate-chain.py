from chainID import sampleChainID
import sys

def main(args):
    import pickle
    import time
    N,M,L = int(args[1]), int(args[2]), int(args[3]) # num vars, var cardinality, no. of limids
    print "generating %d chain-like limids with %d variables taking on %d values..." % (L,N,M),
    sys.stdout.flush()
    start = time.clock()
    for t in range(L):
        #ChanceVars, DecVars, CPT, Strategy, Utility = sampleChainID(N,M)
        exp = sampleChainID(N,M)
        pickle.dump( exp, open( "experiments/exp-%d-%d-%d.p" % (N,M,t), "wb" ) )
    end = time.clock()
    print "done [%gs]" % (end-start)
    print "saved in experiments/exp-%d-%d-L" % (N,M)
        
if __name__ == '__main__':
	import sys
        if len(sys.argv) < 4:
            print "Usage:", sys.argv[0], "num_vars num_var_cardinality num_exps"
            exit(0)
 	main(sys.argv)
