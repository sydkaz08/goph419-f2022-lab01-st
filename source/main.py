import math
import numpy as np 
import matplotlib.pyplot as plt

import Arcsin as AS

# calculate the derivative of equation 17 from lab with respect to ve_v0 and alpha

def partial_derivative(alpha,ve_v0):
    pdtalpha=math.sqrt(1-alpha/(1+alpha)*ve_v0**2)+(1+alpha)*1/2*1/\
        math.sqrt(1-alpha/(1+alpha)*ve_v0**2)*(-ve_v0**2/(1+alpha)**2)
    
    pdtve_v0=(1+alpha)*1/2/math.sqrt(1-alpha/(1+alpha)*ve_v0**2)*\
        (-alpha/(1+alpha)*2*ve_v0)
    return pdtalpha, pdtve_v0

# launch_exact calculates the launch angle for max_alpha and min_alpha
# takes user defined alpha, tol_alpha, and ve_v0
# This computation is done using the built in sin function

def launch_exact(alpha, tol_alpha, ve_v0): #uses the math.asin function to test for error
    min_alpha = (1 - tol_alpha) * alpha
    max_alpha = (1 + tol_alpha) * alpha
    in_rt=1-((min_alpha/(1+min_alpha)*ve_v0**2))
    if in_rt>0:
        try:
            sin_phi=(1+min_alpha)*math.sqrt(in_rt)
            max_angle=math.asin(sin_phi)
            if max_angle>math.pi/2 or max_angle<0:
                max_angle=None
        except:
            max_angle=None
    else:
        max_angle=None
    in_rt=1-((max_alpha/(1+max_alpha)*ve_v0**2))
    if in_rt>0:
        try:
            sin_phi=(1+max_alpha)*math.sqrt(in_rt)
            min_angle=math.asin(sin_phi)
            if min_angle>math.pi/2 or min_angle<0:
                min_angle=None
        except:
            min_angle=None
    else:
        min_angle=None
    return min_angle, max_angle

def launch_approx(alpha, tol_alpha, ve_v0): #using the arcsin() function created using Bowein and Chamberland
    min_alpha = (1 - tol_alpha) * alpha
    max_alpha = (1 + tol_alpha) * alpha
    math_temp=1-((min_alpha/(1+min_alpha)*ve_v0**2))
    if math_temp>0:
        sin_phi=(1+min_alpha)*math.sqrt(math_temp)
        max_angle=AS.arcsin(sin_phi,10)
        if max_angle>math.pi/2 or max_angle<0:
            max_angle=None
    else:
        max_angle=None
    math_temp=1-((max_alpha/(1+max_alpha)*ve_v0**2))
    if math_temp>0:
        sin_phi=(1+max_alpha)*math.sqrt(math_temp)
        min_angle=AS.arcsin(sin_phi,10) #10 chosen to give good approximation
        if min_angle>math.pi/2 or max_angle<0:
            min_angle=None
    else:
        min_angle=None
    return max_angle, min_angle

#To check if the error between the exact and approximate values
#Takes values calculates in launch_angles as x_true and x_approx

def error(la_true,la_approx):
    true_error=abs(la_true-la_approx)
    true_relative_error=true_error/la_true
    print("True Error:",true_error,"\nTrue Relative Error:",true_relative_error)
    
#This function plotts the min and max launch angles over a rangle of alpha values
def fig1():
    ve_v0=2
    tol_alpha=0.04
    num_alpha=int(0.5//0.001+1) #number of alpha 
    alpha=np.linspace(0,0.5,num_alpha)
    min_angle=[0]*num_alpha
    max_angle=[0]*num_alpha
    for i in range(num_alpha):
        max_angle[i],min_angle[i]=launch_approx(alpha[i], tol_alpha, ve_v0)
    plt.plot(alpha,min_angle, label ='Minimum Alpha')
    plt.plot(alpha,max_angle, label ='Maximum Alpha')
    plt.xlabel(chr(945))
    plt.ylabel(chr(966))
    plt.title("Minimum and Maximum Launch Angles with Fluctuating "+str(chr(945))+"\n")
    plt.legend()
    plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab01-st/source/Figures/Figure(1)')
    plt.show()
    
#This function answers part d, plotting the min and max launch angles over a
#range of ve_v0 values
def fig2():
    alpha=0.25
    tol_alpha=0.04
    n_ve_v0=int(2.5//0.001+1)
    ve_v0=np.linspace(0,2.5,n_ve_v0)
    min_angle=[0]*n_ve_v0
    max_angle=[0]*n_ve_v0
    for i in range (n_ve_v0):
        max_angle[i],min_angle[i]=launch_approx(alpha, tol_alpha, ve_v0[i])
    plt.plot(ve_v0,min_angle, label = 'Minimum Alpha')
    plt.plot(ve_v0,max_angle, label ='Maximum Alpha')
    plt.xlabel("ve/v0")
    plt.ylabel(chr(966))
    plt.title("Minimum and Maximum Launch Angles with Fluctuating ve/v0\n")
    plt.legend()
    plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab01-st/source/Figures/Figure(2)')
    plt.show()
    
#This function computes the error in sin_phi, and condition number

#This function computes the error in sin_phi, and condition number
def uncertainty():
    ve_v0=2 
    alpha=0.25
    delta_ve_v0=0.05
    delta_alpha=0.02
    in_rt=1-((alpha/(1+alpha))*(ve_v0**2))
    sin_phi=(1+alpha)*math.sqrt(in_rt)
    pdtalpha,pdtve_v0=partial_derivative(alpha,ve_v0)
    deltasin_phi=abs(pdtalpha)*delta_alpha+abs(pdtve_v0)*delta_ve_v0
    print("\nError of sin(" + str(chr(966)) + "):",str(deltasin_phi))
    

    # From the Notes
    x=math.sqrt(ve_v0**2+alpha**2) #Norm of input values
    f=sin_phi #Function Result
    J=math.sqrt(pdtalpha**2+pdtve_v0**2) #Jacobian Matrix
    CN=x*J/f #solve for CN
    print("\nCondition Number:",CN)

# start() calls previous functions
def start():
    max_angle2,min_angle2=launch_approx(0.25, 0.02, 2)  #APPROX when alpha=0.25, ve_v0=2, alpha_tol=2%
    print("Minumum Angle (APPROX):",str(min_angle2)+"\nMaximum Angle (APPROX):",max_angle2)
    min_angle1,max_angle1=launch_exact(0.25, 0.02, 2)  #TRUE when alpha=0.25, ve_v0=2, alpha_tol=2%
    print("Minumum Angle (EXACT): "+str(min_angle1)+"\nMaximum Angle (EXACT): "+ str(max_angle1)+"\n")
    print("For min_angle")
    error(min_angle1,min_angle2)
    print("\nFor max angle")
    
    error(max_angle1,max_angle2)
    fig1()
    fig2()
    uncertainty()
    
start()
