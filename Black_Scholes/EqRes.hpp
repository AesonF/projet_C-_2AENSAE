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

//classe "portefeuille" ?
class BSparams{
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

#endif /* EqRes_hpp */
