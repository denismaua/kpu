# implements Factor class

from math import exp, log

LOGZERO = -1000.0 #float("-inf")

class Factor:
    """ Class Factor """
    #__slots__ = ['scope','stride','dimension','zeroValue','values','id']
    def __init__(self, variableList, defaultValue=LOGZERO, zeroValue=LOGZERO):
        self.scope = variableList        
        self.stride = {}
        self.dimension = 1
        self.zeroValue = zeroValue
        if len(self.scope) == 1:
            self.stride[self.scope[0]] = 1
            self.dimension = self.scope[0].num_states
        else:
            for variable in self.scope:
                self.stride[variable] = self.dimension
                self.dimension = self.dimension*variable.num_states
        self.values = self.dimension*[defaultValue]
        self.id = Factor.factors
        Factor.factors += 1

    def getStride(self,variable):
        if variable in self.stride:
            return self.stride[variable]
        else:
            return 0
        
    def setValue(self,assignment,value):
        index = 0
        for i,l in enumerate(assignment):
            index += l*self.getStride(self.scope[i])
        self.values[index] = value

    def clamp(self,variable,assignment):
        dimension = self.dimension/variable.num_states
        values = dimension*[0.0]
        j = 0
        for i in range(self.dimension):
            index = int(1.0*i/self.stride[variable]) % variable.cardinality            
            if index == assignment:
                values[j] = self.values[i]
                j += 1
            # if index != assignment:
            #    self.values[i] = self.zeroValue
        del self.stride[variable]
        self.scope.remove(variable)
        self.dimension = 1
        for variable in self.scope:
            self.stride[variable] = self.dimension
            self.dimension = self.dimension*variable.num_states
        self.values = values

    def getValue(self,assignment):
        index = 0
        for i,l in enumerate(assignment):
            index += l*self.getStride(self.scope[i])
        return self.values[index]
    
    def printOut(self):
        print 'Dimension:', self.dimension
        print 'Scope:', [var.label for var in self.scope]
        assignment = len(self.scope)*[0]
        for i in range(self.dimension):
            for j,var in enumerate(self.scope):
                assignment[j] = int(1.0*i/self.stride[var]) % var.num_states
            print i, assignment, self.values[i]
            
Factor.factors = 0
