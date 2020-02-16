# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 18:39:03 2017

@author: Gerard Higgins
"""

import numpy as np
import sys, os
rootDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(rootDir)
os.chdir(rootDir)
from arc2 import *
from wigner import Wigner3j, Wigner6j
import matplotlib.pylab as plt
from scipy import optimize

def linear_function(x,m,c):
    return m*x+c

polarisability_atomic_units=1.648777e-41
h=6.62607004e-34

atom=Strontium()

n1 = 42
l1 = 0
j1 = 0.5

mj1_range=np.arange(0.5,j1+1)

pol1=[]
pol2=[]

###############################################################################
### Find polarisability using getPolarizability command #######################
###############################################################################
calc = StarkMap(atom)
eFieldRange = np.array([0.0,0.1])*1.e2 # V/m ranges
for mj1 in mj1_range:
    calc.defineBasis(n1, l1, j1, mj1, n1-5, n1+5, l1+1, progressOutput=True )        # use low E-field strengths, only need to consider P-states
    minEfield = eFieldRange[0]
    maxEfield = eFieldRange[1]
    calc.diagonalise(np.linspace(minEfield,maxEfield,100))
    pol1.append(calc.getPolarizability( showPlot=False, minStateContribution=0.9))
pol1=np.asarray(pol1)
pol1_au = pol1/polarisability_atomic_units*h/1e4*1e6
print "m1 = ", mj1_range
print "polarisability from getPolarizability command / (MHz cm^2 / V^2)", pol1
print "polarisability from getPolarizability command / (atomic units  )", pol1_au
#print("%.2f MHz cm^2 / V^2" % (pol1))
#print("%.2e atomic units  " % (pol1_au))
print  


###############################################################################
### Find polarisability using dipole matrix elements ##########################
###############################################################################
pol2=[]

n2_range = np.arange(n1-5,n1+5)
if l1==0:
    l2_range=[1]
else:
    l2_range=[l1-1,l1+1]
sum_of_terms=0.
for mj1 in mj1_range:
    pol2_temp=0.
    for n2 in n2_range:
        for l2 in l2_range:
            if l2==0:
                j2_range=[0.5]
            else:
                j2_range=[l2-0.5,l2+0.5]
            for j2 in j2_range:
                mj2_range=np.arange(-j2,j2+1)
                for mj2 in mj2_range:
                    energy_difference = ( atom.getEnergy(n1,l1,j1) - atom.getEnergy(n2,l2,j2)) / 27.2114    # convert eV to Hartree
                    q=mj2-mj1
                    if np.abs(q)<1.1 and np.abs(j2-j1)<1.1:
                        dipole_matrix_element = atom.getDipoleMatrixElement(n1,l1,j1,mj1,n2,l2,j2,mj2,q)
                        radial_matrix_element = atom.getRadialMatrixElement(n1,l1,j1,n2,l2,j2)
                        if j1>0.51:
                            pol2_temp += dipole_matrix_element**2 / energy_difference * (   -2./3.   + (2.*j1+1) * 2.*np.sqrt(10.*j1*(2.*j1-1.)/(3.*(2.*j1+3.)*(j1+1.)*(2.*j1+1.)))*((-1)**(j1+j2+1))*Wigner6j(j1,j2,1.,1.,2.,j1)*((3.*mj1**2-j1*(j1+1.))/(j1*(2.*j1-1.)))   ) 
                        else:
                            pol2_temp += dipole_matrix_element**2 / energy_difference * (   -2./3.   ) 
#                        pol2_temp += dipole_matrix_element**2 / energy_difference 
                        sum_of_terms+=radial_matrix_element**2 / energy_difference

    pol2.append(pol2_temp)
print "here! %.4e" % sum_of_terms
print "here!", atom.getRadialMatrixElement(42,0,0.5,42,1,0.5)
pol2=np.asarray(pol2)
pol2_au=pol2
pol2 = (pol2)*polarisability_atomic_units/h*1e4/1e6
print "m1 = ", mj1_range
print "polarisability calculated from dipole matrix elements / (MHz cm^2 / V^2)", pol2
print "polarisability calculated from dipole matrix elements / (atomic units  )", pol2_au
print

if j1>0.51:
    x_axis=(3.*mj1_range**2-j1*(j1+1.))/(j1*(2.*j1-1.))
    x_axis_fit=np.linspace(min(x_axis),max(x_axis),num=101)
    plt.plot(x_axis,pol1,"*")
    plt.plot(x_axis,pol2,"*")
    popt1,pcov1=optimize.curve_fit(linear_function,x_axis,pol1)
    popt2,pcov2=optimize.curve_fit(linear_function,x_axis,pol2)
    print "alpha_0 %.2f MHz cm^2/V^2 from getPolarizability command" % popt1[1]
    print "alpha_0 %.2f MHz cm^2/V^2 from dipole matrix elements   " % popt2[1]
    print
    print "alpha_0 %.3e atomic units from getPolarizability command" % (popt1[1]/polarisability_atomic_units*h/1e4*1e6)
    print "alpha_0 %.3e atomic units from dipole matrix elements   " % (popt2[1]/polarisability_atomic_units*h/1e4*1e6)
    print
    print "alpha_1 %.2f MHz cm^2/V^2 from getPolarizability command" % popt1[0]
    print "alpha_1 %.2f MHz cm^2/V^2 from dipole matrix elements   " % popt2[0]
    print
    print "alpha_1 %.3e atomic units from getPolarizability command" % (popt1[0]/polarisability_atomic_units*h/1e4*1e6)
    print "alpha_1 %.3e atomic units from dipole matrix elements   " % (popt2[0]/polarisability_atomic_units*h/1e4*1e6)
    print
    plt.plot(x_axis_fit,linear_function(x_axis_fit,popt1[0],popt1[1]),"-")
    plt.plot(x_axis_fit,linear_function(x_axis_fit,popt2[0],popt2[1]),"--")
    print "ratio in alpha_1", popt1[0]/popt2[0]
    plt.xlim([min(x_axis)-0.2,max(x_axis)+0.2])
#    plt.xlabel("3*m^2-J(J+1)/(J(2J-1))")
    plt.show()