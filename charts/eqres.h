#ifndef EQRES_H
#define EQRES_H

#include <stdio.h>
#include <iostream>
#include "matrix.h"

class Matrix; //déclaration anticipée



//========================== DÉRIVÉ QUELCONQUE ========================== //
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



//========================== OPTION VANILLE EU ========================== //
class Vanilla : public BSparams { //call ou put vanille
public:
    Vanilla(float _K, bool _call, int _n, int _m, float _tmax, float _sigma, float _mu);
    std::vector<float> prices(); //vecteur des valeurs à maturité
    float K; //le strike, compris entre 0 et 1 pour simplifier
    bool call;
};



//========================== PORTEFEUILLE D'OVE ========================== //
class Portfolio {
public:
    unsigned int l; //nombre d'options dans le portefeuille
    std::vector<Vanilla> L; //liste des options vanille composant le portefeuille

    Portfolio(); //portefeuille vide
    Portfolio(std::vector<Vanilla> L, unsigned int _l);
    void add(Vanilla V); //ajout d'une option vanille au portefeuille
    std::vector<float> presentValue(unsigned int prec, int method); //valeur du portefeuille au présent
};



//========================== MATRICES DE PASSAGE ========================== //
void matCreate(BSparams &par, Matrix &Mat); //méthode 1 (il faut l'inverser)
Matrix YMatrix(BSparams &par); //méthode 2 (il ne faut pas l'inverser)



//===================== VALORISATIONS (selon S et t) ====================== //
Matrix BSSol(BSparams &par, std::vector<float> &prices, unsigned int method); //valeur d'un dérivé
Matrix BSSol(Vanilla V, unsigned int method); //valeur d'une option vanille européenne
std::vector<float> PFvalue(Portfolio P); //valeur d'un portefeuille



//========================== VALEUR EXACTE D'OVE ========================== //
void singleSim(float tmax, int n, float S0, float sigma, float mu);
float PricingExplicit(float S, float K, float r, float T, float sigma);

std::vector<float> PricingExplicitS(unsigned int prec, Vanilla V);
// (liste des valeurs obtenues avec PricingExplicit)



//========================= EXEMPLES POUR LA DÉMO ========================= //
std::vector<float> priceExamples(unsigned int number, unsigned int gauge);



#endif // EQRES_H
