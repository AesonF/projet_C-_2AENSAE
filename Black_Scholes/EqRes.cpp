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
#include <cmath>
#include <random>


//BSparams DEFINITION DES FONCTIONS DE CLASSE
BSparams::BSparams(int _n, int _m, float _tmax, float _sigma, float _mu){
    n = _n;
    m = _m;
    tmax = _tmax;
    sigma = _sigma;
    mu = _mu;
}

//Matrix DEFINITION DES FONCTIONS DE CLASSE
void Matrix::allocArray()
{
    p = new float*[m];
    for(int i=0; i<m; i++){
        p[i] = new float[n];
    }
}
Matrix::Matrix(int a, int b){
    n = a;
    m = b;
    allocArray();
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            p[j][i] = 0;
        }
    }
}
void Matrix::load(int i, int j, float a){
    p[j][i] = a;
}
float Matrix::iloc(int i, int j) { //get value (i,j)
    return p[j][i];
}
void Matrix::show(){
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            std::cout << p[j][i] << ',';
        }
        std::cout << std::endl;
    }
}

int Matrix::lin(){return n;}
int Matrix::col(){return m;}

Matrix Matrix::operator *(Matrix &A){
    Matrix Res(n,A.col());
    for(int i=0; i<n; i++){
        for(int j=0; j<A.col(); j++){
            float Res_ij = 0; //valeur (i,j) de la matrice resultante de l'operation ("Res")
            for(int k=0; k<m; k++){
                Res_ij += p[k][i]*A.iloc(k,j);
            }
            Res.load(i,j,Res_ij); //insertion de la valeur dans la matrice Res
        }
    }
    return Res;
}

Matrix::~Matrix(){
    for(int j=0; j<m; j++){
        delete [] p[j];
    }
    delete [] p;
}




//ALGORITHMES DE RESOLUTION

/*La valorisation de la dérivée est une fonction du prix du sous-jacent S et
dutemps t. La fonction résolvant l'EDP de Black-Scholes doit donc prendre en
entrée les paramètres de modèle donnes ci-dessus, travaille avec une matrice
de toutes les valeurs pour tous les prix envisageables (discrets ici) et renvoie
un vecteur des valorisations pour un prix donne.*/

//this one creates the transition matrix for going from t to t+1
void matCreate(BSparams &par, std::vector<float> &prices, Matrix &Mat){
    int m = par.m; float fm = (float)(m);
    int n = par.n;
    float dt = par.tmax/float(n);
    float a = 0.0; //variable utilisée pour remplir Mat
    for(int i=0; i<m; i++){ //parcours des lignes
        for(int j=0; j<n; j++){ //parcours des colonnes
            
            if( i-j<-1 or i-j>1 ){a = 0;}
            
            else if(i==j-1){a = 0.5*par.mu*dt-0.5*par.sigma*par.sigma*fm*fm*dt;}
            
            else if(i==j){a = 1 + par.sigma*par.sigma*fm*fm*dt;}
            
            else if(i==j+1){a = -0.5*par.mu*dt-0.5*par.sigma*par.sigma*fm*fm*dt;}
            
            Mat.load(i,j,a);
        }
    }
}

//"matrix[ std::slice( 2, col, row ) ] = pi;" sets third column to pi
//"matrix[ std::slice( 3*row, row, 1 ) ] = e;" sets fourth row to e



/* BSSol = Black-Shcholes solution : prend en entrée une instance de la classe
BSparams, définie en EqRes.hpp, ainsi que le vecteur des prix initiaux du sous-
jacent, sous forme d'un vecteur de taille par.m */

void BSSol(BSparams &par, std::vector<float> &prices){
    std::vector<float> Sol(par.n,0); //0 = valeur par defaut
    Matrix Mat(par.m, par.n);
    matCreate(par, prices, Mat);
}



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
