//  main.cpp
//  HPC_Q3ab
//
//  Created by Desmond Cheung on 17/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//


#include <iostream>
#include <vector>
#include "TriMatrix.h"

using namespace std;


int main(int argc, const char * argv[]) {
    /*----------------------------- Programme description: --------------------------------------------------------------------------------------------------------*/
    
    double L = atof(argv[1]);
    double Nx = atof(argv[2]);
    double T = atof(argv[3]);
    double Nt = atof(argv[4]);
    double alpha = atof(argv[5]);
    double theta = atof(argv[6]);
    
    double dx = L/Nx;
    double dt = T/Nt;
    
    
    /*----------------------------- Other calculated variables ----------------------------------------------------------------------------------------------------*/
    
    double gamma0 = 0;                  //temperature at front end of bar (BC1)
    double gamma1 = 0;                  //temeprature at rear end of bar (BC2)
    double nu = alpha*dt/dx/dx;         //defines the nu constant
    
    /*----------------------------- Defining preliminary vectors used in the programme ----------------------------------------------------------------------------*/
    
    vector<double> x((Nx+1));           //x-co-ordinates
    vector<double> u0((Nx+1));          //initial heat distribution
    vector<double> * u1, * u2;          //vectors declared for matrix multiplication
    u1 = new vector<double>(Nx+1);      //vetors constructed for starting matrix inversion
    u2 = new vector<double>(Nx+1);      //vector defined for receipient of matrix inversion
    
    
    /*----------------------------- Defining the x-co-ordinate space by creating an x-vector ----------------------------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        x[i]=i*dx;                      //increments the x vector according to dx size
        //cout<<x[i]<<endl;             //check that the x-vector makes sense
    }
    
    /*----------------------------- Defining the u0 vector as the initial heat distribution along the bar ---------------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        u0[i]=x[i]*(1-x[i]);            //increments the x vector according to dx size
        //cout<<u0[i]<<endl;            //check that the x-vector makes sense
    }
    
    /*----------------------------- Constructing the starting initial u vector for the heat distribution ----------------------------------------------------------*/
    
    (*u1)[0]=gamma0;
    (*u1)[Nx]=gamma1;
    
    for (int i=1; i<Nx; i++) {
        (*u1)[i] = u0[i]; //ommit the front and end term which was defined above using gamma0 and gamma1
    }
    
    
    /*----------------------------- Generates tri-matrix to be inverted -------------------------------------------------------------------------------------------*/
    TriMatrix LHS(Nx+1, nu, (-1*theta));         //generates the matrix to be inversion
    TriMatrix RHS(Nx+1, nu, (1-theta));          //generates the multiplication matrix
    
    /*----------------------------- Performs time integration by 1st multiplying the LHS of the equation, then taking the inverse of the 2nd equation -------------*/
    
    cout << (*u1)[(Nx+1)/2] << " ";
    for (int t=0; t<Nt-1; t++) {
        u2 = RHS*u1;
        u1 = LHS/u2;
        cout << (*u1)[(Nx+1)/2] << " ";
    }

    
    return 0;
}
