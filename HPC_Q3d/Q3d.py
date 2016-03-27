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
x = 0.5*L;
Nt = 250;   #fix Nt
dtmin = 0.0001;       #sets the lower dt limit  
dtmax = 0.001;      #sets the upper dt limit
dtincre = 0.00005;  #sets dt range

Nxmin = 20;     #sets the starting Nx
pi = math.pi;
dNx = 2;
Nxmax = 30;
Nx = Nxmin;
theta = 0.5;


#end
k=1;
while (Nx < Nxmax + 1):
    dt = dtmin;
    RMS = [];
    dtarray=[];
    j = 1;
    while (dt < dtmax): #calculates RMS error for each dt
        tarray = [];
        Uarray = [];
        for i in range (Nt): #calculates corresponding analytical solution
            t = i*dt;
            tarray.append(dt);
            Ui = 0;
            for n in range (1,100): #calculates fourier series sum
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
        T = dt*Nt;
        dtarray.append(dt);    
        terminalout = subprocess.check_output(['./Q3d',repr(L), repr(Nx), repr(T), repr(Nt), repr(alpha), repr(theta)]);
        Uexp = [float(i) for i in terminalout.split()]; #reads the numerical solution 
        rms = numpy.sqrt(numpy.sum(((Uexp-Uarray)**2))/len(Uexp)); #calculate RMS error associated with specific Nt/Nx step 
        RMS.append(rms);
        dt += dtincre;
        print("RMS points loop number: ", j);
        j += 1;
    #end
    seriesname='Nx = ' + repr(Nx);
    plt.plot(dtarray, RMS, label=seriesname); #plots the RMS vs dt graph for specific Nx case
    Nx += dNx;
    print("Space step loop number: ", k);
    k+=1;
#end
    
plt.title('RMS error Scaling with dt size for the Crank-Nicholson Method');
plt.ylabel('RMS error');
plt.xlabel('dt');
plt.legend(loc='best');