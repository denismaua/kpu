# runs experiments with chain IDs
from Factor import Factor
from Variable import Variable

def Dirichlet(N,M):
    """ Sample a N-dimensional vector from a symmetric Dirichlet with parameter 1/M """
    from random import gammavariate
    samples = [gammavariate(1.0/M,1) for _ in xrange(N)]
    sum_samples = sum(samples)
    return [x/sum_samples for x in samples]

def sampleChainID(N,M):
    """ Sample a chain ID with N chance nodes and variable cardinalities M """
    from random import random
    ChanceVars = [ Variable('C'+str(i),M) for i in range(N) ]
    DecVars = [ Variable('D'+str(i),M) for i in range(N) ]
    CPT = N*[ None ]
    DecCPT = N*[ None ]
    CPT[0] = Factor( [ChanceVars[0],DecVars[0]] ) # P(C[0]|D[0])
    for offset in range(M):            
        # P(C[i]|C[i-1]=0,D[i]=0)
        CPT[0].values[M*offset:M*(offset+1)] = Dirichlet(M,M)
    for i in range(1,N):
        CPT[i] = Factor( [ChanceVars[i],ChanceVars[i-1],DecVars[i]] ) # P(C[i]|C[i-1],D[i])
        offset = 0
        for offset in range(M*M):
            # P(C[i]|C[i-1]=0,D[i]=0)
            CPT[i].values[M*offset:M*(offset+1)] = Dirichlet(M,M)
            offset+=1
    for i in range(N):
        DecCPT[i] = Factor([DecVars[i]], defaultValue=1.0/M)
    U = Factor([ ChanceVars[N-1] ])
    for i in range(M):
        U.values[i] = random()
    return ChanceVars, DecVars, CPT, DecCPT, U

def main(args):
    N,M,NUMEXP = int(args[1]), int(args[2]), int(args[3]) # num vars, var cardinality, no. of experiments to run
    import spu, spu2, spu3, kpu, mpu, sys, akpu
    import time
    from math import log
    
    # sample network
    MPUttime, MPUvalue = NUMEXP*[0.0], NUMEXP*[0.0]
    SPUttime, SPUvalue = NUMEXP*[0.0], NUMEXP*[0.0]
    SPU2ttime, SPU2value = NUMEXP*[0.0], NUMEXP*[0.0]
    SPU3ttime, SPU3value = NUMEXP*[0.0], NUMEXP*[0.0]
    KPUttime, KPUvalue = NUMEXP*[0.0], NUMEXP*[0.0]
    aKPUttime, aKPUvalue = NUMEXP*[0.0], NUMEXP*[0.0]
    aMPUttime, aMPUvalue = NUMEXP*[0.0], NUMEXP*[0.0]
    start = time.clock()
    for t in range(NUMEXP):
        print "### %d ##############################################" % (t+1)
        sys.stdout.flush()
        ChanceVars, DecVars, CPT, Strategy, Utility = sampleChainID(N,M)
        MPU = 0
        if N*log(M) < 22*log(2): # rule of thumb
            print "running MPU...",
            sys.stdout.flush()
            MPU, MPUtime, maxtables = mpu.run(ChanceVars, DecVars, CPT, Strategy, Utility, verbose=False)
            print "done. [\033[91m%gs\033[0m]" % MPUtime
            MPUttime[t] = MPUtime
            MPUvalue[t] = MPU

            #print "#iterations:", SPUit
            print "Solution: %f" % MPU
            print
        ## print "running SPU...",
        ## sys.stdout.flush()
        ## SPUsol = len(DecVars)*[ None ]
        ## for n in range(len(DecVars)):
        ##     SPUsol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
        ## sys.stdout.flush()
        ## SPU, SPUit, SPUtime = spu.run(ChanceVars, DecVars, CPT, SPUsol, Utility, verbose=False)
        ## print "done. [\033[91m%gs\033[0m]" % SPUtime
        ## print "#iterations:", SPUit
        ## print "Solution: %f [%f]" % (SPU,SPU-MPU)
        ## SPUttime += SPUtime
        ## SPUvalue += SPU
        ## print
        print "running SPU2...",
        sys.stdout.flush()
        SPU2sol = len(DecVars)*[ None ]
        for n in range(len(DecVars)):
            SPU2sol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
        SPU2, SPU2it, SPU2time = spu2.run(ChanceVars, DecVars, CPT, SPU2sol, Utility, verbose=False)
        print "done. [\033[91m%gs\033[0m]" % SPU2time
        print "#iterations:", SPU2it
        print "Solution: %f [%f]" % (SPU2,SPU2-MPU)
        SPU2ttime[t] = SPU2time
        SPU2value[t] = SPU2
        print
        print "running SPU2dominance...",
        sys.stdout.flush()
        SPU3sol = len(DecVars)*[ None ]
        for n in range(len(DecVars)):
            SPU3sol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
        SPU3, SPU3it, SPU3time = spu3.run(ChanceVars, DecVars, CPT, SPU3sol, Utility, verbose=False)
        print "done. [\033[91m%gs\033[0m]" % SPU3time
        print "#iterations:", SPU3it
        print "Solution: %f [%f]" % (SPU3,SPU3-MPU)
        SPU3ttime[t] = SPU3time
        SPU3value[t] = SPU3
        print
        print "running KPU...",
        sys.stdout.flush()
        K = 5 # neighborhood size
        MIt = 5*N # num iterations
        KPUsol = SPU3sol
        ## KPUsol = len(DecVars)*[ None ]
        ## for n in range(len(DecVars)):
        ##     KPUsol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
        KPU, KPUit, KPUtime, max_tables = kpu.run(ChanceVars, DecVars, CPT, KPUsol, Utility, K, verbose=False, maxiterations=MIt)
        print "done. [\033[91m%gs\033[0m]" % KPUtime
        print "#iterations:", KPUit
        print "Solution: %f [%f]" % (KPU,KPU-MPU)
        KPUttime[t] = KPUtime + SPU3time
        KPUvalue[t] = KPU
        #print "running SPU2 on KPU solution...",
        sys.stdout.flush()
        ## KPU2, KPU2it, KPU2time = spu2.run(ChanceVars, DecVars, CPT, KPUsol, Utility, verbose=False)
        ## print "done. [\033[91m%gs\033[0m]" % KPU2time
        ## print "#iterations:", SPU2it
        ## print "Solution: %f [%f]" % (KPU2,KPU2-MPU)
        print
        print "running aKPU...",
        sys.stdout.flush()
        K = 5 # neighborhood size
        MIt = 5*N # num iterations
        aKPUsol = SPU2sol
        aKPU, aKPUit, aKPUtime, aKPUsize = akpu.run(ChanceVars, DecVars, CPT, aKPUsol, Utility, K, verbose=False, maxiterations=MIt, Approximate=True, numtables=30000/M)
        print "done. [\033[91m%gs\033[0m]" % aKPUtime
        print "#iterations:", aKPUit
        print "Solution: %f [%f]" % (aKPU,aKPU-MPU)
        aKPUttime[t] = aKPUtime + SPU3time
        aKPUvalue[t] = aKPU
        sys.stdout.flush()
        print
        print "running aMPU...",
        sys.stdout.flush()
        K = N # neighborhood size
        MIt = 1 # num iterations
        #aKPUsol = SPU2sol
        aMPUsol = len(DecVars)*[ None ]
        for n in range(len(DecVars)):
            aMPUsol[n] = Factor([DecVars[n]], defaultValue=1.0/M)
        aMPU, aMPUit, aMPUtime, aMPUsize = akpu.run(ChanceVars, DecVars, CPT, aMPUsol, Utility, K, verbose=False, maxiterations=MIt, Approximate=True, numtables=30000/M)
        
        print "done. [\033[91m%gs\033[0m]" % aMPUtime
        print "#iterations:", aMPUit
        print "Solution: %f [%f]" % (aMPU,aMPU-MPU)
        aMPUttime[t] = aMPUtime
        aMPUvalue[t] = aMPU
        #print "running SPU2 on KPU solution...",
        sys.stdout.flush()
    ttime = time.clock() - start
    print "Finished in \033[91m%gsec\033[0m" % ttime
    print
    print "N=%d, M=%d, T=%d" % (N,M,NUMEXP)
    print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    print "EXPERIMENT | AKPU                                      | SPU                                       | SPUhull                                   | KPU                                       | AMPU"
    print "           |                 TIME                VALUE |                 TIME                VALUE |                 TIME                VALUE |                TIME                 VALUE |                 TIME                VALUE"
    print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    for t in range(NUMEXP):
        print "%10d | %20s %20s | %20s %20s | %20s %20s | %20s %20s | %20s %20s" % (t, str(aKPUttime[t]), str(aKPUvalue[t]), str(SPU2ttime[t]), str(SPU2value[t]), str(SPU3ttime[t]), str(SPU3value[t]), str(KPUttime[t]), str(KPUvalue[t]), str(aMPUttime[t]), str(aMPUvalue[t]))
    print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    print "   AVERAGE | %20s %20s | %20s %20s | %20s %20s | %20s %20s | %20s %20s" % ( str(sum(aKPUttime)/NUMEXP), str(sum(aKPUvalue)/NUMEXP), str(sum(SPU2ttime)/NUMEXP), str(sum(SPU2value)/NUMEXP), str(sum(SPU3ttime)/NUMEXP), str(sum(SPU3value)/NUMEXP), str(sum(KPUttime)/NUMEXP), str(sum(KPUvalue)/NUMEXP), str(sum(aMPUttime)/NUMEXP), str(sum(aMPUvalue)/NUMEXP))
    
if __name__ == '__main__':
	import sys
        if len(sys.argv) < 4:
            print "Usage:", sys.argv[0], "num_vars num_var_cardinality num_exps"
            exit(0)
 	main(sys.argv)
