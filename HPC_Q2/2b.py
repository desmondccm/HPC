# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 03:18:41 2016

@author: 11cheungc3
"""
import subprocess
import matplotlib
import matplotlib.pylab as plt
import numpy
import numpy as np
import math

alpha = 0.1;
L = 1;
T = 5;
x = 0.5*L;
Nt = 1000;
Nx = 20;
dt = T/Nt;

tarray = [];
Uarray = [];
pi = math.pi;

for i in range (Nt):
    t = i*dt;
    tarray.append(t);
    Ui = 0;
    for n in range (1,100): #this forloop is for calculating the fourier series
        Cn = 8*(math.sin(n*pi/2)**2)/pi**3/n**3;
        insidesin = n*pi*x/L;
        insideexp = -alpha*n**2*pi**2*t/L**2;
        term = Cn * math.sin(insidesin)*math.exp(insideexp);
        Ui = term + Ui;
    #end
    Uarray.append(Ui);
#end 
Uarray = numpy.array(Uarray);
tarray = numpy.array(tarray);

terminalout = [];
terminalout = subprocess.check_output(['./Q2',repr(L), repr(Nx), repr(T), repr(Nt), repr(alpha)]);
Uexp = [float(i) for i in terminalout.split()];

plt.figure(1)
plt.subplot(211)
plt.plot(tarray,Uarray,'r', label='Analytical');
plt.title('Analytical vs Numerical Heat Decay at x=2');
plt.legend(loc='best');
plt.ylabel('U');
plt.subplot(212)
plt.plot(tarray,Uexp, label='Numerical');
plt.ylabel('U');
plt.xlabel('t');
plt.legend(loc='best');

plt.figure(2)
plt.title('Analytical vs Numerical Heat Decay at x = 0.5');
plt.plot(tarray,Uarray, label='Analytical');
plt.plot(tarray,Uexp,'r--', label='Numerical');
plt.ylabel('U');
plt.xlabel('t');
plt.legend(loc='best');

