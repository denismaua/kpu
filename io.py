#!/usr/bin/python

from Variable import Variable
from Factor import Factor, LOGZERO

def read_network_from_file(filename, useLog = True, round2int = False, scaleFactor=100, prefix = ""):
    """ Reads network from file in UAI format and returns a pair Variables, Factors """
    #print "Use log:", useLog
    #print "Scale and Round to Integers:", round2int
    #print "Scale Factor:", scaleFactor
    if useLog:
    	from math import log
    infile = open(filename)
    state = 0 # automata initial state
    linecount = 0
    for line in infile:
        linecount += 1
        line = line.lstrip()
        line = line.rstrip()
        if len(line) > 0 and line[0] != '#' and line[0] != '%': # skip blank and commented lines        
            if state == 0:
                # parse header
                if line != "MARKOV" and line != "BAYES":
                    print "In line", linecount
                    print "Expected 'MARKOV' or 'BAYES' file header, found", line
                    exit(1)    
                state += 1
            elif state == 1:
                # parse number of variables
                if not line.isdigit():
                    print "In line", linecount
                    print "Expected number of variables, found", line
                    exit(1)
                num_vars = int(line)
                state += 1
            elif state == 2:
                # parse variable cardinalities
                line = line.split()
                if len(line) != num_vars:
                    print "In line", linecount
                    print "Expected", num_vars, "digits, found", len(line)
                    exit(1)
                    
                Variables = [Variable(prefix + 'X'+str(i+1), int(line[i])) for i in range(num_vars)]
                state += 1
            elif state == 3:
                # parse number of cliques
                if not line.isdigit():
                    print "In line", linecount
                    print "Expected number of cliques, found", line
                    exit(1)
                num_cliques = int(line)
                Factors = []
                state += 1
            elif state == 4:
                # parse factor scopes               
                line = line.split()[1:]
                
                scope = [Variables[int(var)] for var in line]
                scope.reverse()
                f = Factor(scope)
                Factors.append(f)
                num_cliques -= 1
                if num_cliques == 0:
                    state += 1
                    ivalue = 0
                    fdim = 0
            elif state == 5:
                # parse factor values                
                line = line.split()               
                for value in line:
                    if fdim == 0:
                        fdim = int(value)  # factor dimension
                    else:
                        if round2int:
                            if useLog:
                                if float(value) == 0 :
                                    Factors[num_cliques].values[ivalue] = LOGZERO #float("-inf")
                                else:
                                    Factors[num_cliques].values[ivalue] = int(scaleFactor*log(float(value)))
                            else:
                                Factors[num_cliques].values[ivalue] = int(scaleFactor*float(value))  
                        elif useLog:                            
                            if float(value) == 0 :
                                Factors[num_cliques].values[ivalue] = LOGZERO
                            else:
                                Factors[num_cliques].values[ivalue] = log(float(value))
                        else:
    
                            Factors[num_cliques].values[ivalue] = float(value)
                        ivalue += 1
                        if ivalue >= Factors[num_cliques].dimension:
                            ivalue = 0
                            num_cliques += 1
                            fdim = 0
            else:
                print "Error in line", linecount
                exit(1)                           

    return Variables, Factors

def read_assignment_from_file(filename):
    """ Reads full assignments of the variables from file in the UAI MPE format """
    infile = open(filename)
    state = 0 # automata initial state
    num_assigns = 0
    linecount = 0
    assign = []
    for line in infile:
        linecount += 1
        line = line.rstrip()
        if len(line) > 0: # skip blank lines
            if state == 0:
                # parse header
                if line != "MPE":
                    print "In line", linecount
                    print "Expected 'MPE' file header, found", line
                    exit(1)    
                state += 1
            elif state == 1:
                # parse number of assignments
                if not line.isdigit():
                    print "In line", linecount
                    print "Expected number of assignments, found", line
                    exit(1)
                num_assigns = int(line)
                state += 1
            elif state == 2:
                # parse assignments
                line = line.split()
                num_vars = int(line[0])
                line = line[1:]                
                if len(line) != num_vars:
                    print "In line", linecount
                    print "Expected", num_vars, "values, found", len(line)
                    exit(1)
                assign.append([int(i) for i in line])
            else:
                print "Error in line", linecount
                exit(1)                
    return assign

def read_evidence_from_file(filename):
    """ Reads evidence from file in UAI format and returns a dictionary {Variable: evidence} """
    infile = open(filename)
    state = 0
    num_obs = 0
    evidence = {}
    for line in infile:
        line = line.lstrip()
        line = line.rstrip()
        if len(line) > 0 and line[0] != '#' and line[0] != "%":
            if state == 0:
                if not line.isdigit():
                    print "Expected number of observations, found", ' '.join(line)
                    exit(1)
                num_obs = int(line)
                state = num_obs
            else:
                line = line.strip()
                line = line.split()
                for i in range(len(line)/2):
                    evidence[int(line[2*i])] = int(line[2*i+1])
                    state -= 1
    
    return evidence

def save_to_uai_format(filename,Variables,Factors):
    """ Saves a mode to file in UAI format """
    out = open(filename,'w')
    # header
    out.write('MARKOV\n')
    # num vars
    out.write(str(len(Variables))+'\n')
    # var cadinalities
    out.write(' '.join([ str(v.cardinality) for v in Variables ])+'\n')
    # number of cliques
    out.write(str(len(Factors))+'\n')
    # cliques (one per line)
    for f in Factors:
        scope = f.scope[:]
        scope.reverse() # UAI reads from right-to-left, while Factors index configurations from left-to-right
        out.write(str(len(scope))+' '+' '.join([str(Variables.index(v)) for v in scope]) + '\n')                  
    # factor values
    for f in Factors:
        out.write(str(f.dimension)+' '+' '.join([str(v) for v in f.values])+'\n')
    out.close()
