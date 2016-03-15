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

class TriMatrix                     //declare a Tri-Matrix class, a special type of matrix that consist of a non-zero diagonal, upper diagonal and lower diagonal.
{
public:
    int dim;
    double diag, updiag, lowdiag;   //used for matrix generation and matrix properties
    vector<double> vectorin;        //used for matrix multiplication.

    Matrix Matrixgen(void)          //define a function that generates a matrix
    {
        Matrix M(dim, Row(dim));
        for (int i=0; i<dim; i++) {
            M[i][i]=diag;           //filling middle diagonal
            M[i][i+1]=updiag;       //filling upper diagonal
            M[i][i-1]=lowdiag;      //filling lower diagonal
        }
        M[0][0]=1;                  //middle diagonal always starts 1
        M[dim][dim]=1;              //middle diagonal always ends with 1
        M[dim][dim-1]=0;            //lower diagonal always ends with 0
        M[1][2]=0;                  //upper diagonal always starts with 0
        return M;
    }
    
    vector<double> Matrixmult(void)         //defines a function that multiplies matrices to vectors
    {
        Matrix A = Matrixgen();     //Calls the matrix generation function to generate a Matrix
        vector<double> vectorout(dim); //initialize an output vector
        double sum = 0;             //initilizes a sum variable representing each row
        for (int i=0; i<dim; i++) { //defines a forloop that writes member to vector
            for (int j=0; j<dim; j++) {
                double dum=A[i][j]*vectorin[j];
                sum=dum+sum;
            }
            vectorout[i]=sum;
        }
    return vectorout;
    }

};

int main() {
    /*----------------------------- programme description: -------------------------------------------------------------------------------------*/
    cout << "this codes HPC code 1 in MATLAB. The final Vector containing the heat distribution of the bar is displayed" << endl;
    cout << "Timestep is 0.001s, Time domain is 10s, alpha = 1, bar end temperatures = 0 K" << endl;
    
    /*----------------------------- Declare variables used for the problem ---------------------------------------------------------------------*/
    
    double alpha = 1;               //heat conductivity
    double dt = 0.001;              //time step size
    int Nx = 20;                    //space steps
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
    Mat.dim = Nx-1;
    Mat.diag = 1-2*nu;
    Mat.updiag = nu;
    Mat.lowdiag = nu;
    Matrix Matt=Mat.Matrixgen();
    
    /*------------------------ Conduct the actual time-integration process as specified by the matrix multiplying method --------------------------*/
    
    Mat.vectorin = u1;
    for (double t=0; t<T; t+=dt) {
        u2=Mat.Matrixmult();
        Mat.vectorin = u2;
    }
    
    return 0;
}
