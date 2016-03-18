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
    
    cout << "this codes HPC code for question 3c and d. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, Time domain is 10s, alpha = 1, bar end temperatures = 0 K" << endl;
    
    /*----------------------------- Define variables via user input -----------------------------------------------------------------------------------------------*/
    

    double alpha = -1;
    
    
    while (alpha < 0) {                 //heat conductivity, makes sure expression is bigger than 0
        cout << "alpha (suggest 1 > alpha > 0): ";
        cin >> alpha;
    }
    
    cout << endl;
    cout << "dt (suggest smaller than 0.001): ";
    double dt;
    cin >> dt;                          //time step size
    
    cout << endl;
    cout << "Nx (suggest ~20): ";
    double Nx;
    cin >> Nx;                          //Space steps
    
    
    cout << endl;
    cout << "L (suggest 1): ";
    double L;
    cin >> L;                           //bar length
    
    cout << endl;
    cout << "T (suggest O~10s): ";
    double T;
    cin >> T;                           //run time
    
    cout << endl;
    cout << "theta (0<theta<1): ";
    double theta;
    cin >> theta;
    
    /*----------------------------- Other calculated variables ----------------------------------------------------------------------------------------------------*/
    
    double gamma0 = 0;                  //temperature at front end of bar (BC1)
    double gamma1 = 0;                  //temeprature at rear end of bar (BC2)
    double dx = L/Nx;                   //defines space step size
    double nu = alpha*dt/dx/dx;         //defines the nu constant
    double Nt = T/dt;
    
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
        u0[i]=x[i]/(1-x[i]);            //increments the x vector according to dx size
        //cout<<u0[i]<<endl;            //check that the x-vector makes sense
    }
    
    /*----------------------------- Constructing the starting initial u vector for the heat distribution ----------------------------------------------------------*/
    
    (*u1)[0]=gamma0;
    (*u1)[Nx]=gamma1;
    
    for (int i=1; i<Nx; i++) {
        (*u1)[i] = u0[i]; //ommit the front and end term which was defined above using gamma0 and gamma1
    }
    
    
    /*----------------------------- Generates tri-matrix to be inverted -------------------------------------------------------------------------------------------*/
    TriMatrix LHS(Nx+1, nu, (-1*theta));         //generates the matrix inversion
    TriMatrix RHS(Nx+1, nu, (1-theta));          //
    
    /*----------------------------- Performs time integration by 1st multiplying the LHS of the equation, then taking the inverse of the 2nd equation -------------*/
    
    for (int t=0; t<Nt; t++) {
        u2 = RHS*u1;
        u2 = LHS/u2;
        u1 = u2;
    }
    
    cout << endl;
    cout << "This is the resulting heat vector after time: " << T <<"s"<<endl;
    for (int i=0; i<(*u1).size(); i++){
        cout << "x="<<x[i]<<": "<<(*u1)[i] << endl;
    }
    
    return 0;
}
