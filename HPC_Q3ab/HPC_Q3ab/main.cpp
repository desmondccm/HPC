//
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


int main() {
    /*----------------------------- Programme description: --------------------------------------------------------------------------------------------------------*/
    
    cout << "this codes HPC code 1 in MATLAB. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, Time domain is 10s, alpha = 1, bar end temperatures = 0 K" << endl;
    
    /*----------------------------- Define variables via user input -----------------------------------------------------------------------------------------------*/
    
    cout << "alpha: ";
    double alpha;                       //heat conductivity
    cin >> alpha;
    
    cout << endl;
    cout << "dt: ";
    double dt;
    cin >> dt;                          //time step size
    
    cout << endl;
    cout << "Nx: ";
    double Nx;
    cin >> Nx;                          //Space steps
    
    
    cout << endl;
    cout << "L: ";
    double L;
    cin >> L;                           //bar length
    
    cout << endl;
    cout << "T: ";
    double T;
    cin >> T;                           //run time
    
    /*----------------------------- Other calculated variables ----------------------------------------------------------------------------------------------------*/
    
    double gamma0 = 0;                  //temperature at front end of bar (BC1)
    double gamma1 = 0;                  //temeprature at rear end of bar (BC2)
    double dx = L/Nx;                   //defines space step size
    double nu = alpha*dt/dx/dx;         //defines the nu constant
    
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
    
    for (int i=0; i<Nx+1; i++) {        //fill vector based on the size of size of the vector
        u0[i]=x[i]/(1-x[i]);            //increments the x vector according to dx size
        //cout<<u0[i]<<endl;            //check that the x-vector makes sense
    }
    
    /*----------------------------- Constructing the starting initial u vector for the heat distribution ----------------------------------------------------------*/
    
    (*u1)[0]=gamma0;
    (*u1)[Nx]=gamma1;
    
    for (int i=1; i<Nx; i++) {
        (*u1)[i] = u0[i]; //ommit the front and end term which was defined above using gamma0 and gamma1
    }
    
    /*----------------------------- Declare and initialize the 3 vectors that define the TriMatrix according to project brief -------------------------------------*/
    
    vector<double> * updiag;             //declares the upper diagonal via a pointer
    vector<double> * lowdiag;            //declares the main diagonal via a pointer
    vector<double> * maindiag;               //declares the lower diagonal via a pointer
    
    maindiag = new vector<double>(Nx+1);    //constructs a new vector at address diag to store the matrix's diagonal
    lowdiag = new vector<double>(Nx);   //constructs a new vector at address lowdiag to store the matrix's lower diagonal
    updiag = new vector<double>(Nx);    //constructs a new vector at address updiag to store the matrix's upper diagonal
    
    /*----------------------------- Constructs the diagonal vector ------------------------------------------------------------------------------------------------*/
    
    for (int i=1; i<Nx; i++) {
        (*maindiag)[i] = (2 * nu) + 1;      //writes 1+2nu into the middle diagonal except first and last entry
    }
    (*maindiag)[0] = 1;                     //wrties last entry and first entry as 1
    (*maindiag)[Nx] = 1;
    
    /*----------------------------- Constructs the upper diagonal vector ------------------------------------------------------------------------------------------*/
    
    for (int i=1; i<Nx-1; i++) {
        (*updiag)[i] = -1 * nu;             //writes minus nu into the upper diagonal except first and last entry
    }
    (*updiag)[0] = 0;
    (*updiag)[Nx-1] = -1 * nu;
    
    /*----------------------------- Constructs the upper diagonal vector ------------------------------------------------------------------------------------------*/
    
    for (int i=1; i<Nx-1; i++) {
        (*lowdiag)[i] = -1 * nu;            //writes minus nu into the lower diagonal except first and last entry
    }
    (*lowdiag)[Nx-1] = 0;
    (*lowdiag)[0] = -1 * nu;
    
    return 0;
}

