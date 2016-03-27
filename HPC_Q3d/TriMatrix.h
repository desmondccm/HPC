//
//  TriMatrix.h
//  HPC_Q3d
//
//  Created by Desmond Cheung on 15/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//  The trimatrix is defined in such a way that the initiation of the Matrix consists of inputting three vectors specified by the pointers a, b and c, where they
//  correpond to the lower, mid and upper diagonal of the Matrix.

#ifndef CLASS_TriMatrix
#define CLASS_TriMatrix
#include <vector>
using namespace std;


class TriMatrix {
    vector<double> *diag, *uDiag, *lDiag;                               //stores the 3 vectors for the diagonals as pointers to their memory location
    
    
    
    
    
public:                                                                 //Enables public access using the constructor function and other accessor functions
    
    /*--------------------------------------------- Constructor for constructing the trimatrix by defining the diagonal vectors ----------------------------------*/
    
    
    //assumes each entry in the diagonal vectors are unique, this constructor is a universal trimatrix constructor.
    TriMatrix(vector<double> *a, vector<double> *b, vector<double> *c):
    diag(b), uDiag(c), lDiag(a) {};
    
    //special case for this coursework only. This assumes that the diagonal vectors consist of repeating members of specific forms.
    TriMatrix(int dim, double nu, double arg) {
        diag = new vector<double>(dim);
        uDiag = new vector<double>(dim-1);
        lDiag = new vector<double>(dim-1);
        
        for (int i = 0; i < diag->size(); i++) {
            (*diag)[i] = 1-2*nu*arg;
        }
        for (int i = 0; i < uDiag->size(); i++) {
            (*uDiag)[i] = nu*arg;
        }
        for (int i = 0; i < lDiag->size(); i++) {
            (*lDiag)[i] = nu*arg;
        }
        
        (*uDiag)[0] = 0;
        (*lDiag)[dim-2] = 0;
        (*diag)[0] = 1;
        (*diag)[dim-1] = 1;
    }
    
    
    
    /*--------------------------------------------- Matrix multiplication function accessable by public ----------------------------------------------------------*/
    
    vector<double> *operator *(vector<double> *U){                      //initiates a matrix mult function that multiplies matrices with vectors by the * operator
        int dim = diag->size();
        vector<double> *u2;
        vector<double> dum1(dim), dum2(dim), dum3(dim);
        u2 = new vector<double>(dim);
        
        for (int i = 0; i < dim; i++) {                                 //for adding the middle diagonal's contribution to the output vector
            dum1[i] = (*U)[i] * (*diag)[i];
        }
        
        for (int i = 1; i < dim ; i++) {                                //for adding the lower diagonal's contribution to the output vector
            dum2[i] = (*U)[i-1] * (*lDiag)[i-1];
        }
        
        for (int i = 0; i < (dim - 1); i++) {                           //for adding the upper diagonal's contribution to the output vector
            dum3[i] = (*U)[i+1] * (*uDiag)[i];
        }
        
        for (int i = 0; i < (dim - 1); i++) {
            (*u2)[i] = dum1[i] + dum2[i] + dum3[i];
        }
        (*u2)[0] = dum1[0] + dum3[0];
        (*u2)[dim-1] = dum1[dim-1] + dum2[dim-1];
        
        delete U;
        
        return u2;                                                      //returns multiplcation product
    }
    
    
    /*--------------------------------------------- Matrix printing function accessable by public ---------------------------------------------------------------*/
    
    void display(){
        cout << "the matrix consists of: " << endl;
        
        cout << "Upper diagonal: " << endl;                             //printing the upper diagonal
        for (int i=0; i < (*uDiag).size(); i++) {
            cout << (*uDiag)[i] << "  ";
        }
        cout << endl;
        cout << "Main diagonal: " << endl;                              //printing main diagonal
        for (int i=0; i < (*diag).size(); i++) {
            cout << (*diag)[i] << "  ";
        }
        cout << endl;
        cout << "Lower diagonal: " << endl;                             //printing main diagonal
        for (int i=0; i < (*lDiag).size(); i++) {
            cout << (*lDiag)[i] << "  ";
        }
    }
    
    /*--------------------------------------------- Matrix inversing function accessable by public ---------------------------------------------------------------*/
    
    vector<double> *operator/(vector<double> *y){                       //uses the forward slash operator to perform inversing function
        int dim = diag->size();
        vector<double> *x, diagPrime(dim), yPrime(dim), newLdiag(dim);
        x = new vector<double>(dim);
        
        for (int i = 1; i < dim; i++) {
            newLdiag[i] = (*lDiag)[i-1];
        }
        
        yPrime[0]=(*y)[0];
        diagPrime[0] = (*diag)[0];
        
        for (int i = 1; i < dim ; i++) {                                 //performs forward substitution phase
            double m = newLdiag[i]/diagPrime[i-1];
            diagPrime[i] = (*diag)[i] - m*(*uDiag)[i-1];
            yPrime[i] = (*y)[i] - m*(yPrime)[i-1];
        }
        
        
        
        (*x)[dim-1] = yPrime[dim-1]/diagPrime[dim-1];                   //performs backwards substitution phase
        for (int i=(dim-2); i>=0; i--){
            (*x)[i] = (yPrime[i] - (*uDiag)[i] * (*x)[i+1])/diagPrime[i];
        }
        
        delete y;
        return x;                                                       //returns outputs vector
    }
    
    ~TriMatrix(){};
    
};

#endif
