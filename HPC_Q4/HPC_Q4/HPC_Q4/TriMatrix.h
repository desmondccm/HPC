//
//  TriMatrix.h
//  HPC
//
//  Created by Desmond Cheung on 15/03/2016.
//  Copyright (c) 2016 cmc213. All rights reserved.
//  The trimatrix is defined in such a way that the initiation of the Matrix consists of inputting three vectors specified by the pointers a, b and c, where they
//  correpond to the lower, mid and upper diagonal of the Matrix.

#ifndef CLASS_TriMatrix
#define CLASS_TriMatrix
#include <vector>
#include <Accelerate/Accelerate.h>
//#include <cblas.h>
//#include <lapacke.h>

using namespace std;

class TriMatrix {
    int dim;
    double *diag, *uDiag, *lDiag;                               //stores the 3 vectors for the diagonals as pointers to their memory location
    double *matarray;                                           //stores the matrix as a vector for BLAS usage

    
    
    
    
public:                                                                 //Enables public access using the constructor function and other accessor functions
    
    /*--------------------------------------------- Constructor for constructing the trimatrix by defining the diagonal vectors ----------------------------------*/
    
    
    //assumes each entry in the diagonal vectors are unique, this constructor is a universal trimatrix constructor.
    TriMatrix(double *a, double *b, double *c):
    diag(b), uDiag(c), lDiag(a) {};
    
    //special case for this coursework only. This assumes that the diagonal vectors consist of repeating members of specific forms.
    TriMatrix(int length, double nu, double arg) {
        diag = new double[length];
        uDiag = new double[length-1];
        lDiag = new double[length-1];
        
        for (int i = 0; i < length; i++) {
            diag[i] = 1-2*nu*arg;
        }

        for (int i = 0; i < length - 1; i++) {
            uDiag[i] = nu*arg;
        }
        for (int i = 0; i < length - 1; i++) {
            lDiag[i] = nu*arg;
        }
        
        uDiag[0] = 0;
        lDiag[length-2] = 0;
        diag[0] = 1;
        diag[length-1] = 1;
        dim = length;
    }
    
    
    
    /*--------------------------------------------- Matrix multiplication function accessable by public ----------------------------------------------------------*/
    
    
    /*--------------------------------------------- Matrix printing function accessable by public ---------------------------------------------------------------*/
    
    void display(){
        cout << "the matrix consists of: " << endl;
        
        cout << "Upper diagonal: " << endl;                             //printing the upper diagonal
        for (int i=0; i < dim - 1; i++) {
            cout << uDiag[i] << "  ";
        }
        cout << endl;
        cout << "Main diagonal: " << endl;                              //printing main diagonal
        for (int i=0; i < dim; i++) {
            cout << diag[i] << "  ";
        }
        cout << endl;
        cout << "Lower diagonal: " << endl;                             //printing main diagonal
        for (int i=0; i < dim - 1; i++) {
            cout << lDiag[i] << "  ";
        }
    }
    
    
    /*---------------------------------------------- This converts the matrix to a vector -----------------------------------------------------------------------*/
    void Mat2vec(){
        double matarray[dim * dim];
        matarray[0] = diag[0];                                         //creating the first 2 entries for the converted vector
        matarray[1] = lDiag[0];
        int j = 0;
        
        for (int i = 2; i < dim; i++){                         //creating the next (Nx-1) entries for the matrix
            matarray[i] = 0;
        }
        for (int i = dim; i < dim * dim - dim; i = i + dim){                         //creating a list of vectors
            matarray[i + j] = uDiag[i/dim-1];
            matarray[i+1 + j] = diag[i/dim];
            matarray[i+2 + j] = uDiag[i/dim];
            j++;
        }
        matarray[dim*dim-1] = diag[dim-1];
        matarray[dim*dim-2] = uDiag[dim-2];
        matarray[dim*dim-dim-1] = 0;

    }
    
    /*---------------------------------------------- This multiplies the matrix using the BLAS routine -------------------------------------------------------------*/
    
    void multiblas(double *x, int length){
        double y[length];                                //define the x-array as an array

        double alpha, beta;
        alpha = 1;
        beta = 0;
    
        cblas_dgemv(CblasColMajor, CblasNoTrans, length, length, alpha, matarray, length, x, 1, beta, y, 1); //BLAS HERE!!!!!!!!!!!!!!!!!

        x = y;
        
    }
    
    /*---------------------------------------------- This inverse the matrix using LAPACK dgtsv function ---------------------------------------------------------------*/
    void inlapack(double *x, int length){

        int nrhs = 1;
        int info;
        
        dgtsv_(&length, &nrhs, lDiag, diag, uDiag, x, &length, &info);
        
        if (info != 0){
            cout << "DGTSV routine error! Error code: " << info <<endl;
        }
    }
    
        
    
};

#endif
