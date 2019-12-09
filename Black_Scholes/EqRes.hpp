//
//  EqRes.hpp
//  Black_Scholes
//
//  Created by Yacine Falaki and Aeson Feehan.
//
//  DESCRIPTION:
//  algorithms for the resolution of the Black-Scholes equation

#ifndef EqRes_hpp
#define EqRes_hpp

#include <stdio.h>
#include <iostream>
#include "Matrix.hpp"

class Matrix; //déclaration anticipée

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

void matCreate(BSparams &par, std::vector<float> &prices, Matrix &Mat);

void BSSol(BSparams &par, std::vector<float> &prices);

void singleSim(float tmax, int n, float S0, float sigma, float mu);

#endif /* EqRes_hpp */
