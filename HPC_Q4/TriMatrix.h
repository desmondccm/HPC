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
    double *diag, *uDiag, *lDiag, *matarray;                  //stores the diagonals, the matrix arrays
    double *dl, *d, *du, *du2;                                //stores the factorized diagonals and the super diagonal
    int *ipiv;
    


    
    
    
    
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
        matarray = new double[dim * dim];
        matarray[0] = diag[0];                                          //creating the first 2 entries for the converted vector
        matarray[1] = lDiag[0];
        int j = 0;
        
        for (int i = 2; i < dim; i++){                                  //creating the next (Nx-1) entries for the matrix
            matarray[i] = 0;
        }
        for (int i = dim; i < dim * dim - dim; i = i + dim){            //creating a list of vectors
            matarray[i + j] = uDiag[i/dim-1];
            matarray[i+1 + j] = diag[i/dim];
            matarray[i+2 + j] = uDiag[i/dim];
            j++;
        }
        matarray[dim*dim-1] = diag[dim-1];
        matarray[dim*dim-2] = uDiag[dim-2];
        matarray[dim*dim-dim-1] = 0;
    }
    
    void disparray(){
        for (int i = 0; i < dim * dim; i++){
            cout << matarray[i] << "  ";
        }
    }
    /*---------------------------------------------- This multiplies the matrix using the BLAS routine -------------------------------------------------------------*/

    double *multiblas(double *x){
        double *y = new double[dim];                                    //define the x-array as an array
        
        double alpha, beta;
        alpha = 1;
        beta = 0;
        
        cblas_dgemv(CblasColMajor, CblasNoTrans, dim, dim, alpha, matarray, dim, x, 1, beta, y, 1); //BLAS HERE!!!!!!!!!!!!!!!!!
        
        //I want to assign the value of y to x.
        delete x;
        return y;
    }
    
    /*---------------------------------------------- This inverse the matrix using LAPACK dgtsv function ---------------------------------------------------------------*/
    void matfactor(){
        int info;
        dl = new double[dim-1];
        d = new double[dim];
        du = new double[dim-1];
        du2 = new double[dim-2];
        ipiv = new int[dim];
        
        copy(lDiag, lDiag + dim-1, dl);
        copy(diag, diag + dim, d);
        copy(uDiag, uDiag + dim-1, du);
        
    
        dgttrf_(&dim, dl, d, du, du2, ipiv, &info);
        
        if (info != 0){
            cout << "dgttrd routine (factoring) error! Error code: " << info <<endl;
        }
    }
    
    double *matsolve(double *b){
        int nrhs = 1;
        int info;
        char trans = 'N';
        
        dgttrs_(&trans, &dim, &nrhs, dl, d, du, du2, ipiv, b, &dim, &info);
        
        if (info != 0){
            cout << "DGTTRS routine (solve) error! Error code: " << info <<endl;
        }
        
        return b;
    }
        
    
};

#endif
