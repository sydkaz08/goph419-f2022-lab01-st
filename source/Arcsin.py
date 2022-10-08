import math 
import numpy as np 
import matplotlib.pyplot as plt



#This function is for employing the Browein and Chamberland equations for arcsin
#Takes input value x, N is for the level of precision
def arcsin(x,N):
    arcsin_sqrt=0 #to not get unboundlocal error
    for n in range(N): #N is for level of precision
        n=n+1
        approx_1=(2*x)**(2*n)
        approx_2=n**2*(math.factorial(2*n)/((math.factorial(n))**2))
        arcsin_sqrt=arcsin_sqrt+0.5*approx_1/approx_2
    arcsin=math.sqrt(arcsin_sqrt)
    return arcsin