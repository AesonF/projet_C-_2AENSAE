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
}

int Matrix::lin(){return n;}
int Matrix::col(){return m;}

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

Matrix::~Matrix(){
    for(int j=0; j<m; j++){
        delete [] p[j];
    }
    delete [] p;
}
