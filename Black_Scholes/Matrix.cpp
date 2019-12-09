//
//  Matrix.cpp
//  Black_Scholes
//
//  Created by Yacine Falaki and Aeson Feehan.
//
//  DESCRIPTION:
//  implementation of a matrix class for solving the Black-Scholes equation
//  with a finite differences-type method

#include "Matrix.hpp"
#include <iostream>
#include <vector>

void Matrix::allocArray()
{
    p = new float*[m];
    for(int i=0; i<m; i++){
        p[i] = new float[n];
    }
}

Matrix::Matrix(int a, int b){
    n = a;
    m = b;
    allocArray();
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            p[j][i] = 0;
        }
    }
}

Matrix::Matrix(int a, int b, std::vector<float> &A){
    n = a;
    m = b;
    allocArray();
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            p[j][i] = A[i*m + j];
        }
    }
}

void Matrix::load(int i, int j, float a){
    p[j][i] = a;
}

float Matrix::load(int i, int j) { //get value (i,j)
    return p[j][i];
}

void Matrix::show(){
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            std::cout << p[j][i] << ',';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int Matrix::lin(){return m;}

int Matrix::col(){return n;}

Matrix Matrix::operator *(Matrix &A){
    Matrix Res(n,A.col());
    for(int i=0; i<n; i++){
        for(int j=0; j<A.col(); j++){
            float Res_ij = 0; //valeur (i,j) de la matrice resultante de l'operation ("Res")
            for(int k=0; k<m; k++){
                Res_ij += p[k][i]*A.load(k,j);
            }
            Res.load(i,j,Res_ij); //insertion de la valeur dans la matrice Res
        }
    }
    return Res;
}
void Matrix::operator =(Matrix &A){
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            p[j][i] = A.load(i,j);
        }
    }
}
void Matrix::operator *=(Matrix &A){
    Matrix Res(n,A.col());
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            float Res_ij = 0; //valeur (i,j) de la matrice resultante de l'operation ("Res")
            for(int k=0; k<m; k++){
                Res_ij += p[k][i]*A.load(k,j);
            }
            Res.load(i,j,Res_ij); //insertion de la valeur dans la matrice Res
        }
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            p[j][i] = Res.load(i,j);
        }
    }
}

Matrix::~Matrix(){
    for(int j=0; j<m; j++){
        delete [] p[j];
    }
    delete [] p;
}

Matrix Matrix::copy(){
    Matrix Cp(n,m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            Cp.load(i,j,p[j][i]);
        }
    }
    return Cp;
}

Matrix id(int a){
    Matrix M(a,a);
    for(int i=0; i<a; i++){
        M.load(i,i,1);
    }
    return M;
}

Matrix quickExp(Matrix &A, int a){
    int dim = A.col();
    Matrix C = id(dim);
    Matrix B = A.copy();
    Matrix temp(dim,dim); //ici, n=m, sinon l'exponentiation n'a pas de sens.
    int b = 1;
    while(a>b){
        while(a%2==0 and a>0){
            Matrix temp = B*B;
            B = temp;
            b *= 2;
            a -= b;
        }
        C *= B;
        a -= b;
    }
    return C;
}
