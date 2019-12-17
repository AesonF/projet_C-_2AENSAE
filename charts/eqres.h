#ifndef EQRES_H
#define EQRES_H

#include <stdio.h>
#include <iostream>
#include "matrix.h"

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

void matCreate(BSparams &par, Matrix &Mat);

Matrix BSSol(BSparams &par, std::vector<float> &prices);

void singleSim(float tmax, int n, float S0, float sigma, float mu);

std::vector<float> priceExample1(unsigned int gauge);
std::vector<float> priceExample2(unsigned int gauge);


#endif // EQRES_H
