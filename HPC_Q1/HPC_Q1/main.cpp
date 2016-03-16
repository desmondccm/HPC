//
//  main.cpp
//  HPC_Q1
//
//  Created by Desmond Cheung on 14/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//

#include <iostream>
#include <vector>
#include "TriMatrix.h"

using namespace std;


int main() {
    /*----------------------------- Programme description: --------------------------------------------------------------------------------------------------------*/
    
    cout << "this codes HPC code 1 in MATLAB. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, Time domain is 10s, alpha = 1, bar end temperatures = 0 K" << endl;
    
    /*----------------------------- Declare variables used for the problem ----------------------------------------------------------------------------------------*/
    
    double alpha = 1;                   //heat conductivity
    double dt = 0.001;                  //time step size
    double Nx = 20;                     //space steps
    double L = 1;                       //bar length
    double gamma0 = 0;                  //temperature at front end of bar (BC1)
    double gamma1 = 0;                  //temeprature at rear end of bar (BC2)
    double T = 10;                      //run time
    double dx = L/Nx;                   //defines space step size
    double nu = alpha*dt/dx/dx;         //defines the nu constant
    
    /*----------------------------- Defining preliminary vectors used in the programme ----------------------------------------------------------------------------*/
    
    vector<double> x((Nx+1));           //x-co-ordinates
    vector<double> u0((Nx+1));          //initial heat distribution
    vector<double> u1((Nx+1));          //vector defined for matrix multiplication
    vector<double> u2((Nx+1));          //vector defined for receipient of matrix multiplcation
    

    /*----------------------------- Defining the x-co-ordinate space by creating an x-vector ----------------------------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        x[i]=i*dx;                      //increments the x vector according to dx size
        //cout<<x[i]<<endl;             //check that the x-vector makes sense
    }
    
    /*----------------------------- Defining the u0 vector as the initial heat distribution along the bar ---------------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        u0[i]=x[i]/(1-x[i]);            //increments the x vector according to dx size
        //cout<<u0[i]<<endl;            //check that the x-vector makes sense
    }
    
    /*----------------------------- Constructing the starting initial u vector for the heat distribution ----------------------------------------------------------*/
   
    u1[1]=gamma0;                       //first entry of u1 is 0
    u1[(Nx+1)]=gamma1;                  //last entry of u1 is also 0
    
    for (int i=1; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        u1[i]=u0[i];                    //increments the x vector according to dx size
        //cout<<u1[i-1]<<endl;          //check that the x-vector makes sense
    }
                                        //cout<<u1[Nx+1]<<endl;
    
    /*----------------------------- Declare and initialize the 3 vectors that define the TriMatrix according to project brief -------------------------------------*/
    
    vector<double> *updiag;             //declares the upper diagonal's via a pointer
    vector<double> *lowdiag;            //declares the upper diagonal's via a pointer
    vector<double> *diag;               //declares the upper diagonal's via a pointer
    
    diag = new vector<double>(Nx+1);    //constructs a new vector at address diag to store the matrix's diagonal
    lowdiag = new vector<double>(Nx+1); //constructs a new vector at address lowdiag to store the matrix's lower diagonal
    updiag = new vector<double>(Nx+1);  //constructs a new vector at address updiag to store the matrix's upper diagonal
    
    /*----------------------------- Constructs the diagonal vector ------------------------------------------------------------------------------------------------*/

    for (int i=1; i<Nx; i++) {
        (*diag)[i] = 1 - (2 * nu);      //writes 1-2nu into the middle diagonal except first and last entry
    }
    (*diag)[0] = 1;                     //wrties last entry
    (*diag)[Nx] = 1;
    
    /*----------------------------- Constructs the upper diagonal vector ------------------------------------------------------------------------------------------*/
    
    for (int i=1; i<Nx-1; i++) {
        (*updiag)[i] = nu;              //writes nu into the upper diagonal except first and last entry
    }
    (*updiag)[0] = 0;
    
    /*----------------------------- Constructs the upper diagonal vector ------------------------------------------------------------------------------------------*/
    
    for (int i=1; i<Nx-1; i++) {
        (*updiag)[i] = nu;              //writes nu into the lower diagonal except first and last entry
    }
    (*updiag)[Nx-1] = 0;
    
    /*----------------------------- Conduct the actual time-integration process as specified by the matrix multiplying method -------------------------------------*/
    
    MMM = new TriMatrix(updiag, lowdiag, diag);

    return 0;
}
