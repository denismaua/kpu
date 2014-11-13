#!/bin/python

import pickle
import sys
from Factor import Factor
from Variable import Variable

def main(filename):
    t = 0
    ChanceVars, DecVars, CPT, DecCPT, U = pickle.load( open(filename, "rb") )
    N = len(ChanceVars)
    M = ChanceVars[0].cardinality
    print "LIMID\n" + str(len(ChanceVars)) + " " + str(len(DecVars)) + " 1"
    v_id = {}
    i = 0
    for var in ChanceVars:
        v_id[var.id] = i
        print str(var.cardinality),
        i += 1
    print
    for var in DecVars:
        v_id[var.id] = i            
        print str(var.cardinality),
        i += 1
    print
    for j,var in enumerate(ChanceVars):
        print str(len(CPT[j].scope)-1),
        for pa in CPT[j].scope[1:]:
            print str(v_id[pa.id]),
        print
    for var in DecVars:
        print "0"
    print "1 " + str(v_id[U.scope[0].id])

    for i,var in enumerate(ChanceVars):
        print str(len(CPT[i].values)) 
        print '',
        for val in CPT[i].values:
            print str(val),
        print
    print str(len(U.values))
    print '',
    for val in U.values:
        print str(val),
    print   

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:", sys.argv[0], "filename"
        exit(0)
    main(sys.argv[1])
