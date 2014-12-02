# this is example 1 in the paper, which generates a case with 2^n non-dominated potentials
n = 3

from SetFactor import SetFactor, SetFactorProduct, SetFactorSumProduct, SetFactorSumBProduct
from Variable import Variable

A = Variable("A",2)
B = Variable("B",2)
C = Variable("C",2)
D = [ Variable("D"+str(i),2) for i in range(n) ]
PA = SetFactor([A]+D)
PB = SetFactor([B]+D)
PA.addEmptyTable()
PB.addEmptyTable()
for i in range(2**n):
    c = [(i/2**j)%2 for j in range(n)]
    #print i, c, sum(float(c[j])/2**(j+1) for j in range(n))
    x = sum(float(c[j])/2**(j+1) for j in range(n))
    PA.setValue(0,[1]+c,x)
    PA.setValue(0,[0]+c,1-x)
    x = sum(float(1-c[j])/2**(j+1) for j in range(n))
    PB.setValue(0,[1]+c,x)
    PB.setValue(0,[0]+c,1-x)
PC = SetFactor([C,A,B])
PC.addEmptyTable()
PC.setValue(0,[1,0,0],0.0)
PC.setValue(0,[1,1,0],0.5)
PC.setValue(0,[1,0,1],0.5)
PC.setValue(0,[1,1,1],1.0)
PC.setValue(0,[0,0,0],1.0)
PC.setValue(0,[0,1,0],0.5)
PC.setValue(0,[0,0,1],0.5)
PC.setValue(0,[0,1,1],0.0)
PD = [ SetFactor([Di]) for Di in D ]
for P in PD:
    P.addTable([0,1])
    P.addTable([1,0])
U = SetFactor([C])
U.addTable([0.0, float(2**(n+1))])
P0 = SetFactorSumBProduct(PD+[PA,PB],D)
P1 = SetFactorSumProduct(U,PC,[C])
P3 = SetFactorProduct(P1,P0)
P3.printOut() 

