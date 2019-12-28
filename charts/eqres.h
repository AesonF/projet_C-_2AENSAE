#ifndef EQRES_H
#define EQRES_H

#include <stdio.h>
#include <iostream>
#include "matrix.h"

class Matrix; //déclaration anticipée

//classe des paramètres de modèle
class BSparams {
public:
    BSparams(int _n, int _m, float _tmax, float _sigma, float _mu);
    BSparams();
    int n;
    int m;
    float tmax;
    float sigma;
    float mu;
};

class Vanilla : public BSparams { //call ou put vanille
private:
    float K; //le strike, compris entre 0 et 1 pour simplifier
    bool call;
public:
    Vanilla(float _K, bool _call, int _n, int _m, float _tmax, float _sigma, float _mu);
    std::vector<float> prices(unsigned int prec);
};

class Portfolio {
private:
    unsigned int l; //nombre d'options dans le portefeuille
    std::vector<Vanilla> L; //liste des options vanille composant le portefeuille
public:
    Portfolio(); //portefeuille vide
    Portfolio(std::vector<Vanilla> L, unsigned int _l);
    void add(Vanilla V); //ajout d'une option vanille au portefeuille
    std::vector<float> presentValue(unsigned int prec, int method); //valeur du portefeuille au présent
};

void matCreate(BSparams &par, Matrix &Mat);

Matrix YMatrix(BSparams &par); //Yacine's matrix

Matrix BSSol(BSparams &par, std::vector<float> &prices, unsigned int method);

void singleSim(float tmax, int n, float S0, float sigma, float mu);

std::vector<float> priceExamples(unsigned int number, unsigned int gauge);


#endif // EQRES_H
