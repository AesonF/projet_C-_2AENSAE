//
//  EqRes.cpp
//  Black_Scholes
//
//  Created by Aeson Feehan on 30/11/2019.
//  Copyright © 2019 Aeson Feehan. All rights reserved.
//

#include "EqRes.hpp"
#include <valarray> //type matriciel (ca ne marche pas)
#include <vector>



/* BSSol = Black-Shcholes solution : prend en entrée une instance de la classe
BSparams, définie en EqRes.hpp, ainsi que le vecteur des prix initiaux du sous-
jacent, sous forme d'un vecteur de taille par.m
 
La valorisation de la dérivée est une fonction du prix du sous-jacent S et
dutemps t. La fonction résolvant l'EDP de Black-Scholes doit donc prendre en
entrée les paramètres de modèle donnes ci-dessus, travaille avec une matrice
de toutes les valeurs pour tous les prix envisageables (discrets ici) et renvoie
un vecteur des valorisations pour un prix donne.*/

//this one creates the transition matrix for going from t to t+1
void matCreate(BSparams &par, std::vector<float> &prices, std::valarray<float> &Mat){
    int m = par.m; float fm = (float)(m);
    int n = par.n;
    float dt = par.tmax/float(n);
    float a = 0; //variable utilisée pour remplir Mat
    for(int i=0; i<m; i++){ //parcours des lignes
        for(int j=0; j<n; j++){ //parcours des colonnes
            
            if( i-j<-1 or i-j>1 ){a = 0;}
            
            else if(i==j-1){a = 0.5*par.mu*dt*-0.5*par.sigma*par.sigma*fm*fm*dt;}
            
            else if(i==j){a = 1+ par.sigma*par.sigma*fm*fm*dt;}
            
            else if(i==j+1){a = -0.5*par.mu*dt*-0.5*par.sigma*par.sigma*fm*fm*dt;}
            
            Mat[i,j] = a;
        }
    }
}

void dispMat(std::valarray<float> &Mat){
    
}

void BSSol(BSparams &par, std::vector<float> &prices){
    std::vector<float> Sol(par.n,0); //0 = valeur par defaut
    std::valarray<float> Mat(par.m, par.n);
    matCreate(par, prices, Mat);
    //dispMat(&Mat);
}

//"matrix[ std::slice( 2, col, row ) ] = pi;" sets third column to pi
//"matrix[ std::slice( 3*row, row, 1 ) ] = e;" sets fourth row to e
