//
//  main.cpp
//  HPC_Q1
//
//  Created by Desmond Cheung on 14/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//

#include <iostream>
#include <math.h>
using namespace std;

int main() {
    cout << "this codes HPC code 1 in MATLAB. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, run time is 10s, alpha = 1, bar end temperatures = 0 K and ";
    
    //Declare variables used for the problem
    double alpha = 1; //heat conductivity
    double dt = 0.001; //time step size
    int Nx = 20; //space steps
    double L = 1; //bar length
    double gamma0 = 0; //temperature at front end of bar (BC1)
    double gamma1 = 0; //temeprature at rear end of bar (BC2)
    double T = 10; //run time
    double Nt = T/dt; //number of time steps
    double dx = L/Nx;
    double nu = alpha*dt/dx/dx;
    
    
    
    return 0;
}
