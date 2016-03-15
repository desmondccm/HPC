//
//  main.cpp
//  HPC_Q1
//
//  Created by Desmond Cheung on 14/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//

#include <iostream>
#include <vector>

using namespace std;
typedef vector<double> Row;         //declare a single row
typedef vector<Row> Matrix;         //Define a matrix as a colum of rows

using namespace std;
typedef vector<double> Row;         //declare a single row
typedef vector<Row> Matrix;         //Define a matrix as a colum of rows


class Vector: public std::vector<double> {
public:
    Vector(int size): std::vector<double>(size) {}
    void fill(double val) {
        for(int i=0; i<this->size(); i++)
            this[i] = val;
    }
};


class TriMatrix {
    Vector *diag, *suDiag, *sbDiag;
public:
    TriMatrix(int Nx, double nu) {
        int Nx2 = Nx * Nx;
        diag = new Vector(Nx);
        suDiag = new Vector(Nx2-1);
        sbDiag = new Vector(Nx2-1);
        diag->fill(1-2*nu);
    }
    vector<double> Matrixmult(vector<double> U){
        vector<double> u2;
        
        
        return u2;
    }
};


int main() {
    /*----------------------------- programme description: -------------------------------------------------------------------------------------*/
    cout << "this codes HPC code 1 in MATLAB. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, Time domain is 10s, alpha = 1, bar end temperatures = 0 K" << endl;
    
    /*----------------------------- Declare variables used for the problem ---------------------------------------------------------------------*/
    
    double alpha = 1;               //heat conductivity
    double dt = 0.001;              //time step size
    double Nx = 20;                 //space steps
    double L = 1;                   //bar length
    double gamma0 = 0;              //temperature at front end of bar (BC1)
    double gamma1 = 0;              //temeprature at rear end of bar (BC2)
    double T = 10;                  //run time
    double dx = L/Nx;               //defines space step size
    double nu = alpha*dt/dx/dx;     //defines the nu constant
    
    /*----------------------------- Defining all vectors and matrices used in the programme -----------------------------------------------------*/
    
    vector<double> x((Nx+1));       //x-co-ordinates
    vector<double> u0((Nx+1));      //initial heat distribution
    vector<double> u1((Nx+1));      //vector defined for matrix multiplication
    vector<double> u2((Nx+1));      //vector defined for receipient of matrix multiplcation
    

    /*----------------------------- Defining the x-co-ordinate space by creating an x-vector -----------------------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {    //fill vector based on the size of size of the vector
        x[i]=i*dx;    //increments the x vector according to dx size
        //cout<<x[i]<<endl; //check that the x-vector makes sense
    }
    
    /*----------------------------- Defining the u0 vector as the initial heat distribution along the bar ----------------------------------------*/
    
    for (int i=0; i<Nx+1; i++) {    //fill vector based on the size of size of the vector
        u0[i]=x[i]/(1-x[i]);        //increments the x vector according to dx size
        //cout<<u0[i]<<endl;        //check that the x-vector makes sense
    }
    
    /*----------------------------- constructing the starting initial u vector for the heat distribution ------------------------------------------*/
    u1[1]=gamma0;
    u1[(Nx+1)]=gamma1;
    
    for (int i=1; i<Nx+1; i++) {    //fill vector based on the size of size of the vector
        u1[i]=u0[i];                //increments the x vector according to dx size
                                    //cout<<u1[i-1]<<endl; //check that the x-vector makes sense
    }
                                    //cout<<u1[Nx+1]<<endl;
    
    /*------------------------ Construct the Tri-diagonal Matrix Mat by defining it as an object within class TriMatrix ---------------------------*/
    TriMatrix Mat;
    Mat.dim = (Nx+1);
    Mat.diag = (1-2*nu);
    Mat.updiag = nu;
    Mat.lowdiag = nu;
    Matrix Matt=Mat.Matrixgen();
    
    /*------------------------ Conduct the actual time-integration process as specified by the matrix multiplying method --------------------------*/
    
    Mat.vectorin = u1;
    for (double t=0; t<T; t+=dt) {
        u2=Mat.Matrixmult();
        Mat.vectorin = u2;
        cout<<u2[1];
    }

    return 0;
}
