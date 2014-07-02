#!/usr/bin/python

class Variable:
    """ Class Variable """
    #__slots__ = ['label', 'num_states', 'cardinality', 'id']
    def __init__(self, name, num_states):
        self.label = name
        self.num_states = num_states
        self.cardinality = num_states
        self.id = Variable.variables
        Variable.variables += 1

    def __as_immutable__(self):
        return self.id

Variable.variables = 0

