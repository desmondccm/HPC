//  main.cpp
//  HPC_Q4
//
//  Created by Desmond Cheung on 17/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//


#include <iostream>
#include <math.h>
#include <vector>
#include "TriMatrix.h"
#include <Accelerate/Accelerate.h> //this header is Mac OS specific because the code is developed in XCode, which comes with Accelerate framework from apple. Please disregard this header and compile using the linux header to compile the programme.
//#include <cblas.h>
//#include <lapacke.h>

using namespace std;


int main() {
    /*----------------------------- Programme description: --------------------------------------------------------------------------------------------------------*/
    
    cout << "this codes HPC code for question 3c and d. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, Time domain is 10s, alpha = 1, bar end temperatures = 0 K" << endl;
    
    /*----------------------------- Define variables via user input -----------------------------------------------------------------------------------------------*/
    
    
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
    
    
    cout << endl;
    cout << "theta (0<theta<1): ";
    double theta;
    cin >> theta;
    
    /*----------------------------- Other calculated variables ----------------------------------------------------------------------------------------------------*/
    
    double gamma0 = 0;                  //temperature at front end of bar (BC1)
    double gamma1 = 0;                  //temeprature at rear end of bar (BC2)
    double nu = alpha*dt/dx/dx;         //defines the nu constant

    cout << "The CFL factor (nu) is: " << nu << endl;
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
    RHS.display();
    
    /*----------------------------- Performs time integration by 1st multiplying the LHS of the equation, then taking the inverse of the 2nd equation -------------*/

    
    cout << endl;

    
    vector<double> a((Nx+1)*(Nx+1));
    
    RHS.Mat2vec();
    
    for (int t = 0; t < Nt; t++){
        u2 = RHS.multiblas(u1, (Nx+1));
        u2 = LHS.inlapack(u2, (Nx+1));
        u1 = u2;
    }


    cout << endl;
    cout << "This is the resulting heat vector after time: " << T <<"s"<<endl;
    for (int i=0; i<(*u1).size(); i++){
        cout << "x="<<x[i]<<": "<<(*u1)[i] << endl;
    }
    
    return 0;
}
