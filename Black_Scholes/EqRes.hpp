//
//  EqRes.hpp
//  Black_Scholes
//
//  Created by Aeson Feehan on 30/11/2019.
//  Copyright © 2019 Aeson Feehan. All rights reserved.
//
// EqRes contains algorithms for the resolution of the Black-Scholes equation.

#ifndef EqRes_hpp
#define EqRes_hpp

#include <stdio.h>
#include <iostream>

//classe des paramètres de modele
class BSparams {
public:
    BSparams(int _n, int _m, float _tmax, float _sigma, float _mu);
    int n;
    int m;
    float tmax;
    float sigma;
    float mu;
};

class Matrix {
private:
    int n;
    int m;
    float **p;
    void allocArray();
public:
    Matrix(int a, int b); //constructeur d'une matrice de zeros
    void load(int i, int j, float a); //insertion d'une valeur a en (i,j)
    float iloc(int i, int j); //valeur a l'emplacement (i,j)
    int lin(); //nombre de lignes
    int col(); //nombre de colonnes
    void show(); //affichage d'une matrice
    Matrix operator *(Matrix &A);
    ~Matrix();
};

void matCreate(BSparams &par, std::vector<float> &prices, Matrix &Mat);

void BSSol(BSparams &par, std::vector<float> &prices);

void singleSim(float tmax, int n, float S0, float sigma, float mu);

#endif /* EqRes_hpp */
