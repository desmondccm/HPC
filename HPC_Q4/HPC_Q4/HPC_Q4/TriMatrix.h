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

using namespace std;

class TriMatrix {
    vector<double> *diag, *uDiag, *lDiag;                               //stores the 3 vectors for the diagonals as pointers to their memory location
    vector<double> *matarray;                                           //stores the matrix as a vector for BLAS usage
    
    
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
        
        delete y;                                                       //cleans up
        return x;                                                       //returns outputs vector
    }
    
    /*---------------------------------------------- This converts the matrix to a vector -----------------------------------------------------------------------*/
    void Mat2vec(){
        vector<double> dum(diag->size()*diag->size());
        dum[0] = (*diag)[0];                                         //creating the first 2 entries for the converted vector
        dum[1] = (*lDiag)[0];
        int j = 0;
        
        for (int i = 2; i < diag->size(); i++){                         //creating the next (Nx-1) entries for the matrix
            dum[i] = 0;
        }
        for (int i = diag->size(); i < (diag->size())*(diag->size())-(diag->size()); i = i + diag->size()){                         //creating a list of vectors
            dum[i + j] = (*uDiag)[i/(diag->size())-1];
            dum[i+1 + j] = (*diag)[i/(diag->size())];
            dum[i+2 + j] = (*uDiag)[i/(diag->size())];
            j++;
        }
        dum[diag->size()*diag->size()-1] = (*diag)[diag->size()-1];
        dum[diag->size()*diag->size()-2] = (*uDiag)[uDiag->size()-1];
        dum[diag->size()*diag->size()-diag->size()-1] = 0;
        
        matarray = new vector<double>(diag->size()*diag->size());

        for (int i = 0; i < (diag->size()*diag->size()); i++){
            (*matarray)[i] = dum[i];
        }

    }
    
    /*---------------------------------------------- This multiplies the matrix using the BLAS routine -------------------------------------------------------------*/
    
    vector<double> *multiblas(vector<double> *xold, int length){
        double x[length];                                //define the diagonal length as the leading dimension of the matrix
        double y[length];                                //define the x-array as an array
        double a[length*length];                         //
        vector<double> *output;
        output = new vector<double>(length);
        
        for (int i = 0; i < length * length; i++){
            a[i]=(*matarray)[i];
            //cout << a[i] << "  ";
        }
        //cout << endl;
        
        copy(xold->begin(), xold->end(), x);
        
        /*for (int i = 0; i < length; i++){
            cout<<x[i]<<endl;
        }
        */
        
        double alpha, beta;
        alpha = 1;
        beta = 0;
    
        cblas_dgemv(CblasColMajor, CblasNoTrans, length, length, alpha, a, length, x, 1, beta, y, 1); //BLAS HERE!!!!!!!!!!!!!!!!!

        
        for (int i = 0; i < length * length; i++){
            (*output)[i]=y[i];
        }
        
        return output;
    }
    
    /*---------------------------------------------- This inverse the matrix using LAPACK dgtsv function ---------------------------------------------------------------*/
    vector<double> *inlapack(vector<double> *xold, int length){
        double x[length];                                //define the diagonal length as the leading dimension of the matrix
        double low[length-1], center[length], up[length-1];
        vector<double> *output;
        output = new vector<double>(length);
        int nrhs = 1;
        int info;
       
        for (int i = 0; i < length - 1; i++){
            low[i] = (*lDiag)[i];
            up[i] = (*diag)[i];
        }
        
        for (int i = 0; i < length; i++){
            center[i]=(*diag)[i];
        }
        
        //cout << endl;
        
        copy(xold->begin(), xold->end(), x);
        
        dgtsv_(&length, &nrhs, low, center, up, x, &length, &info);
        
        for (int i = 0; i < length * length; i++){
            (*output)[i]=x[i];
        }
        
        
        
        
        return output;
        
    }
    
};

#endif
