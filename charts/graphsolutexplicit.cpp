#include <iostream>
#includen"math.h"
#include "matrix.h"
#include "eqres.h"

class Matrix;
class BSparams;

using namespace std;

double Pricingexplicit(double S,double K,double r,double T,double sigma){
    double d1 = ( log(S/K) + (r + sigma*sigma /2)*T)/(sigma*sqrt(T) );
    double d2 = d1 - sigma*sqrt(T);
    double Ke = K*exp(-r*T);
    double N1 = 0.5*erfc(-d1/sqrt(2));
    double N2 = 0.5*erfc(-d2/sqrt(2));
    double C = S*N1 - Ke*N2;
    return C;
}


Matrix PricingexplicitS(Matrix S,double K,double r,double T,double sigma){
Matrix C(par.m,1);
for(int i=0,i<,i++){
    C.load(i,1,Pricingexplicit(S.load(i,1),K , r, T,sigma));

return C;
}

int main()
{
    double n = Pricingexplicit(3850,4100,0.0125,1,0.0168);
}
