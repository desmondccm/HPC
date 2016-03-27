//
//  TriMatrix.h
//  HPC_Q1
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
    
    TriMatrix(vector<double> *a, vector<double> *b, vector<double> *c): //list initialization
    diag(b), uDiag(c), lDiag(a) {};
    
    
    
    /*--------------------------------------------- Matrix multiplication function accessable by public ----------------------------------------------------------*/
    
    vector<double> *operator *(vector<double> *U){                      //initiates a matrix multiplication function * that overloads the traditional * function 
        int dim = diag->size();
        vector<double> *u2;
        vector<double> dum1(dim), dum2(dim), dum3(dim);
        u2 = new vector<double>(dim);
        
        
        for (int i = 0; i < (dim - 1); i++) {                           //for adding the upper diagonal's contribution to the output vector
            dum3[i] = (*U)[i+1] * (*uDiag)[i];
        }

        
        for (int i = 0; i < dim; i++) {                                 //for adding the middle diagonal's contribution to the output vector
            dum1[i] = (*U)[i] * (*diag)[i];
        }
        
        for (int i = 1; i < dim ; i++) {                                //for adding the lower diagonal's contribution to the output vector
            dum2[i] = (*U)[i-1] * (*lDiag)[i-1];
        }
        
        for (int i = 0; i < (dim - 1); i++) {
            (*u2)[i] = dum1[i] + dum2[i] + dum3[i];
        }
        (*u2)[0] = dum1[0] + dum3[0];
        (*u2)[dim-1] = dum1[dim-1] + dum2[dim-1];
        
        delete U;                                                       //frees up memory by freeing up pointer space
        
        return u2;
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
    
    
};

#endif
