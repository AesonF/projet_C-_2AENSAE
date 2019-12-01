//
//  Matrix.hpp
//  Black_Scholes
//
//  Created by Yacine Falaki and Aeson Feehan.
//
//  DESCRIPTION:
//  implementation of a matrix class for solving the Black-Scholes equation
//  with a finite differences-type method

#ifndef Matrix_hpp
#define Matrix_hpp

#include <stdio.h>

class Matrix {
private:
    int n;
    int m;
    float **p;
    void allocArray();
public:
    Matrix(int a, int b); //constructeur d'une matrice de zeros
    void load(int i, int j, float a); //insertion d'une valeur a en (i,j)
    float load(int i, int j); //observation de la valeur en (i,j)
    int lin(); //nombre de lignes
    int col(); //nombre de colonnes
    void show(); //affichage d'une matrice
    Matrix operator *(Matrix &A);
    ~Matrix();
};

#endif /* Matrix_hpp */
