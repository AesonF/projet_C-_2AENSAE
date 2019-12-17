//
//  main.cpp
//  Black_Scholes
//
//  Created by Yacine Falaki and Aeson Feehan.

#include <iostream>
#include <vector>
#include "EqRes.hpp"
#include "Matrix.hpp"

void BSmatTest(){ //test de résolution du modèle de Black-Scholes (pas encore fonctionnel)
    BSparams par(5,5,1,0.1,0.1);
    std::vector<float> prices = priceExample1(4);
    BSSol(par, prices);
}

void matTest(){
    int vints[]={1,1,1,0,1,1,0,0,1};
    std::vector<float> vect(vints, vints + sizeof(vints)/sizeof(int));
    Matrix A(3,3,vect);
    Matrix At = transpose(A);
    Matrix B = A*At;
    B.show();
    Matrix C = inverse(B);
    C.show();
}

void bigMatTest(){
    BSparams par(5,5,1,0.1,0.1);
    Matrix A(5,5);
    matCreate(par,A);
    Matrix B = inverse(A);
    Matrix C = A*B;
    C.show();
    //LaTeXShow(C);
}

int main(int argc, const char * argv[]) {
    //BSmatTest();
    //matTest();
    bigMatTest();
    return 0;
}
