//
//  EqRes.cpp
//  Black_Scholes
//
//  Created by Yacine Falaki and Aeson Feehan.
//
//  DESCRIPTION:
//  algorithms for the resolution of the Black-Scholes equation

#include "EqRes.hpp"
#include "Matrix.hpp"
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
    Matrix Sol(par.m,1,prices);
    /* "On ne peut pas faire des mathématiques en confondant les
        vecteurs et les matrices colonnes."
                              — Sandie Souchet
    */
    Matrix Mat(par.m, par.n);
    matCreate(par, prices, Mat);
    Matrix m = Mat.copy();
    for(int t=0; t<par.tmax; t++){
        m *= Mat;
    }
    Matrix temp = Sol*m;
    Sol = temp;
    Sol.show();
}



/* priceExample1 :
Générateur d'un exemple de vecteur de prix initiaux d'un sous-jacent, pour mettre
en entrée de la fonction BSSol. "gauge" correspond à la finesse du découpage de
l'intervalle [1,2], dans lequel applique la fonction f:x->1+x^2 pour 0<x<1.*/
std::vector<float> priceExample1(int gauge){
    std::vector<float> prices(gauge);
    for(int k=0; k<gauge; k++){
        prices[k] = 1+(float(k)/float(gauge))*(float(k)/float(gauge));
        std::cout << prices[k] << ", ";
    }
    std::cout << std::endl;
    return prices;
}



/* SINGLE SIMULATION :
 Cette première valorisation d'un call européen vanille dans le cadre du
modèle se fait de manière itérative, dans une seule fonction. Arguments :
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
