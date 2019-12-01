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
    int **p;
    
    void allocArray()
    {
        p = new int*[m];
        for(int i=0; i<m; i++){
            p[i] = new int[n];
        }
    }
public:
    Matrix(int a, int b){
        n = a;
        m = b;
        allocArray();
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                p[i][j] = 0;
            }
        }
    }
    
    void load(int i, int j, type a){
        p[i][j] = a;
    }
    
    type iloc(int i, int j) { //get value (i,j)
        return p[i][j];
    }
    
    void show(){
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                std::cout << p[i][j] << ',';
            }
            std::cout << std::endl;
        }
    }
    
    ~Matrix(){
        for(int i = 0; i < m; i++){
            delete [] p[i];
        }
        delete [] p;
    }
};

#endif /* EqRes_hpp */
