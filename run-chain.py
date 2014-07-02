import sys

def main(algo,filenames):
    import pickle
    import time
    from Factor import Factor
    #print "MODEL"+(35*' '), "ALGORITHM            NUMVARS    NUMSTATES  ITERATIONS      NUMTABLES       TIME                 VALUE               "
    for filename in filenames:
        #print "loading model from file '%s'..." % filename,
        print filename.ljust(40), algo.ljust(20),
        sys.stdout.flush()
        start = time.clock()
        ChanceVars, DecVars, CPT, Strategy, Utility = pickle.load( open(filename, "rb"))
        end = time.clock()
        N, M = len(DecVars), DecVars[0].cardinality
        print str(N).ljust(10), str(M).ljust(10),
        sys.stdout.flush()
        #print "done [%gs]" % (end-start)
        sys.stdout.flush()
        if algo == "SPU":
            import spu2
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/DecVars[n].cardinality)
            value, it, runtime = spu2.run(ChanceVars, DecVars, CPT, sol, Utility, verbose=False)
            print str(it).ljust(15),'1'.ljust(15),str(runtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        elif algo == "SPUhull":
            import spu3
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/DecVars[n].cardinality)
            value, it, runtime = spu3.run(ChanceVars, DecVars, CPT, sol, Utility, verbose=False)
            print str(it).ljust(15),'1'.ljust(15),str(runtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        elif algo == "MPU":
            import mpu
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/DecVars[n].cardinality)
            value, runtime, numtables = mpu.run(ChanceVars, DecVars, CPT, sol, Utility, verbose=False)            
            print '1'.ljust(15),str(numtables).ljust(15),str(runtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        elif algo[0].isdigit() and algo[-2:] == "PU":
            import spu3, kpu
            i=0
            while algo[i].isdigit():
                i+=1
            K = int(algo[0:i]) # neighborhood size
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/DecVars[n].cardinality)
            value, it, runtime = spu3.run(ChanceVars, DecVars, CPT, sol, Utility, verbose=False)
            MIt = 10*N # num iterations
            value, it, KPUtime, numtables = kpu.run(ChanceVars, DecVars, CPT, sol, Utility, K, verbose=False, maxiterations=MIt)
            print str(it).ljust(15),str(numtables).ljust(15),str(runtime+KPUtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        elif algo[0] == "a" and algo[1].isdigit() and algo[-2:] == "PU":
            import spu3, akpu
            i=1
            while algo[i].isdigit():
                i+=1
            K = int(algo[1:i]) # neighborhood size
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/DecVars[n].cardinality)
            value, it, runtime = spu3.run(ChanceVars, DecVars, CPT, sol, Utility, verbose=False)
            MIt = 10*N # num iterations
            value, it, KPUtime, numtables = akpu.run(ChanceVars, DecVars, CPT, sol, Utility, K, verbose=False, maxiterations=MIt, Approximate=True, numtables=2**18/M/M)
            print str(it).ljust(15),str(numtables).ljust(15),str(runtime+KPUtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        elif algo == "APU":
            import akpu
            K = N # neighborhood size
            MIt = 1 # num iterations
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
            #value, it, runtime, numtables = akpu.run(ChanceVars, DecVars, CPT, sol, Utility, K, verbose=False, maxiterations=MIt, Approximate=True, numtables=2**20/M/M) # numtables is set so that set-potentials size is at most 2**20
            value, it, runtime, numtables = akpu.run(ChanceVars, DecVars, CPT, sol, Utility, K, verbose=False, maxiterations=MIt, Approximate=True, numtables=5000) # numtables is set so that set-potentials size is at most 2**20
            print str(it).ljust(15),str(numtables).ljust(15),str(runtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        elif algo == "MIX":
            import spu3, akpu
            prevvalue = 0.0
            value = 1.0
            MIt = 100
            runtime = 0.0
            maxtables = 0
            epoch = 0
            sol = len(DecVars)*[ None ]
            for n in range(len(DecVars)):
                sol[n] = Factor([DecVars[n]], defaultValue=1.0/DecVars[n].cardinality)
            while prevvalue != value and epoch <= MIt:
                # epoch cycle
                prevvalue = value
                # run spu until convergence
                value, it, sputime = spu3.run(ChanceVars, DecVars, CPT, sol, Utility, verbose=False)
                runtime += sputime
                K = 10 # neighborhood size
                # run N k-local searches
                #value, it, kputime, numtables = kpu.run(ChanceVars, DecVars, CPT, sol, Utility, K, verbose=False, maxiterations=N)
                value, it, kputime, numtables = akpu.run(ChanceVars, DecVars, CPT, sol, Utility, K, verbose=False, maxiterations=N, Approximate=True, numtables=2**18/M/M) 
                maxtables = max(maxtables,numtables)
                runtime += kputime
                epoch += 1
            print str(epoch).ljust(15),str(maxtables).ljust(15),str(runtime).ljust(20),str(value).ljust(20)
            sys.stdout.flush()
        else:
            print "algorithm %s not available." % algo
            exit(1)
        


if __name__ == '__main__':
	import sys
        if len(sys.argv) < 2:
            print "Usage:", sys.argv[0], "algorithm filename1 [filename2 ...] "
            print "\t algorithm is one of:"
            print "\t\t SPU     (no pruning)"
            print "\t\t SPUhull (with dominance pruning)"
            print "\t\t K-PU    (K is a number specifying the neighborhood size)"
            print "\t\t K-APU   (K is a number specifying the neighborhood size)"
            print "\t\t APU"
            print "\t\t MPU"
            exit(0)
 	main(sys.argv[1],sys.argv[2:])
