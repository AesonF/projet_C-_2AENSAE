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

//classe "portefeuille" ?
class BSparams {
public:
    BSparams(int _n, int _m, float _tmax, float _sigma, float _mu){
        n = _n;
        m = _m;
        tmax = _tmax;
        sigma = _sigma;
        mu = _mu;
    }
    int n;
    int m;
    float tmax;
    float sigma;
    float mu;
};

template <typename type>
class Matrix {
private:
    int n;
    int m;
    type val[]; //valeurs, rangees dans une seule tres longue ligne
public:
    Matrix(int a, int b){
        n = a;
        m = b;
        type V[n*m];
        for(int k=0; k<n*m; k++){
            V[k] = 0;
        }
        val = V;
    }
    
    type iloc(int i, int j) { //get value (i,j)
        return val[i*m + j];
    }
    
    void show(){
        for(int k=0; k<n*m-1; k++){
            std::cout << val[k] << ',';
            if(k%n == 0){
                std::cout << std::endl;
            }
        }
        std::cout << val[n*m-1] << std::endl;
    }
};

#endif /* EqRes_hpp */
