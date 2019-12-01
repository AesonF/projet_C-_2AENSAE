//
//  main.cpp
//  Black_Scholes
//
//  Created by Yacine Falaki and Aeson Feehan.

#include <iostream>
#include <vector>
#include "EqRes.hpp"

void BSmatTest(){
    BSparams par(5,5,1,0.01,0.02);
    std::vector<float> prices{0.2,0.4,0.6,0.8,1.0};
    BSSol(par, prices);
}

void matTest(){
    Matrix m1(2,2);
    m1.load(0,0,1); m1.load(1,0,0); m1.load(0,1,0); m1.load(1,1,1);
    m1.show();
    Matrix m2(2,1);
    m2.load(0,0,2); m2.load(1,0,1);
    Matrix m3 = m1*m2;
    m3.show();
}

int main(int argc, const char * argv[]) {
    BSmatTest();
    matTest();
    return 0;
}
