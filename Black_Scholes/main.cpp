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
    BSparams par(100,100,1,0.1,0.1);
    std::vector<float> prices = priceExample1(100);
    BSSol(par, prices);
}

void matTest(){
    std::vector<float> vect{1,0,0,1,1,0,1,1,1};
    Matrix A(3,3,vect);
    Matrix B = quickExp(A,4);
    B.show();
}

int main(int argc, const char * argv[]) {
    BSmatTest();
    //matTest();
    return 0;
}
