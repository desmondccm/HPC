//
//  test.cpp
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
    
    cout << "HPC Question 1: The final Vector containing the heat distribution of the bar is displayed" << endl;
    
    /*----------------------------- Declare variables used for the problem ----------------------------------------------------------------------------------------*/
    
    double alpha = -0.1;
    
    
    while (alpha < 0) {                 //heat conductivity, makes sure expression is bigger than 0
        cout << "alpha (suggest 1 > alpha > 0): ";
        cin >> alpha;
    }
    
    
    
    cout << endl;
    cout << "L (suggest 1): ";
    double L = 1;
    cin >> L;                           //bar length
    
    
    cout << endl;
    cout << "Nx (suggest ~20): ";
    double Nx;
    cin >> Nx;                          //Space steps
    
    double dx = L/Nx;
    
    double dtlim = 0.5 * dx * dx / alpha;      //computes the maximum time step using CFL criteria
    
    cout << endl;
    cout << "T (suggest O~10s): ";
    double T;
    cin >> T;                           //run time
    
    int Ntlim = T/dtlim + 1;
    
    cout << "By computing the maximum CFL criteria, the suggested maximum number of timestep should be: " << Ntlim << endl;
    
    cout << endl;
    cout << "number of time steps: ";
    double Nt;
    cin >> Nt;                          //time step size
    
    double dt = T/Nt;
    
    if (Nt < Ntlim){
        cout << "Entered time step is less than recommended time step and slution may not converge!!" <<endl;
    }
    

    
    /*----------------------------- Other calculated variables ----------------------------------------------------------------------------------------------------*/
    
    double gamma0 = 0;                  //temperature at front end of bar (BC1)
    double gamma1 = 0;                  //temeprature at rear end of bar (BC2)
    double nu = alpha*dt/dx/dx;         //defines the nu constant
    
    cout << "The CFL factor (nu) is: " << nu << endl;
    
    /*----------------------------- Defining preliminary vectors used in the programme ----------------------------------------------------------------------------*/
    
    vector<double> x((Nx+1));           //x-co-ordinates
    vector<double> u0((Nx+1));          //initial heat distribution
    vector<double> * u1, * u2;          //vectors declared for matrix multiplication
    u1 = new vector<double>(Nx+1);      //vetors constructed for starting matrix multiplication
    u2 = new vector<double>(Nx+1);      //vector defined for receipient of matrix multiplcation
    
    
    /*----------------------------- Defining the x-co-ordinate space by creating an x-vector ----------------------------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        x[i]=i*dx;                      //increments the x vector according to dx size
        //cout<<x[i]<<endl;             //check that the x-vector makes sense
    }
    
    /*----------------------------- Defining the u0 vector as the initial heat distribution along the bar ---------------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {
        u0[i]=x[i]*(1-x[i]);
    }
    
    /*----------------------------- Constructing the starting initial u vector for the heat distribution ----------------------------------------------------------*/
    
    (*u1)[0]=gamma0;
    (*u1)[Nx]=gamma1;
    
    for (int i=1; i<Nx; i++) (*u1)[i] = u0[i]; //ommit the front and end term which was defined above using gamma0 and gamma1
    
    /*----------------------------- Declare and initialize the 3 vectors that define the TriMatrix according to project brief -------------------------------------*/
    
    vector<double> * updiag;             //declares the upper diagonal via a pointer
    vector<double> * lowdiag;            //declares the main diagonal via a pointer
    vector<double> * diag;               //declares the lower diagonal via a pointer
    
    diag = new vector<double>(Nx+1);    //constructs a new vector at address diag to store the matrix's diagonal
    lowdiag = new vector<double>(Nx);   //constructs a new vector at address lowdiag to store the matrix's lower diagonal
    updiag = new vector<double>(Nx);    //constructs a new vector at address updiag to store the matrix's upper diagonal
    
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
    (*updiag)[Nx-1] = nu;
    
    /*----------------------------- Constructs the upper diagonal vector ------------------------------------------------------------------------------------------*/
    
    for (int i=1; i<Nx-1; i++) {
        (*lowdiag)[i] = nu;              //writes nu into the lower diagonal except first and last entry
    }
    (*lowdiag)[Nx-1] = 0;
    (*lowdiag)[0] = nu;
    
    /*----------------------------- Conduct the actual time-integration process as specified by the matrix multiplying method -------------------------------------*/
    
    TriMatrix triMatrix(lowdiag, diag, updiag);
    
    triMatrix.display();
    
    for(double i=0; i<(T-dt); i+=dt) {
        u2 = triMatrix*u1;
        u1 = u2;
    }
    
    cout << endl;
    cout << "This is the resulting heat vector after time: " << T <<"s"<<endl;
    for (int i=0; i<(*u1).size(); i++){
        cout << "x="<<x[i]<<": "<<(*u1)[i] << endl;
    }
    return 0;
}
