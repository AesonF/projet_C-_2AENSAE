//
//  main.cpp
//  Black_Scholes
//
//  Created by Aeson Feehan on 30/10/2019.
//  Copyright © 2019 Aeson Feehan. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <random>

/* SINGLE SIMULATION :
 Cette première simulation du cours d'un produit dérivé dans le cadre du modèle de
 Black-Scholes se fait de manière itérative, dans une seule fonction. Arguments :
 - tmax  = temps de la simulation
 - n     = nombre de subdivisions temporelles
 - S0    = prix initial du sous-jacent
 - sigma = volatilité du sous-jacent
 - mu    = tendance du sous-jacent
*/
void singleSim(float tmax, int n, float S0, float sigma, float mu)
{
    float T[n]; T[0] = 0 ;  //Axe des abscisses (temps)
    float S[n]; S[0] = S0; //Axe des ordonnées (prix du sous-jacent)
    float S1 = 0;
    float dt = tmax/n;
    
    std::default_random_engine generator;
    std::normal_distribution<float> distribution(0.0,dt);
    
    float norm0 = 0;
    float norm1 = 0;
    for(int i=1; i<n+1 ; i++){
        norm0 = norm1;
        norm1 = distribution(generator);
        S1 = S0*(mu*dt + sigma*(norm1+norm0) + 1);
        S[i] = S1;
        T[i] = tmax*i/n;
        S0 = S1;
    }
}


int main(int argc, const char * argv[]) {
    singleSim(10,1000,100,0.3,0.1);
}
