#include "eqres.h"
#include "matrix.h"

#include <vector>
#include <cmath>
#include <random>

class Matrix;

//========================== DÉRIVÉ QUELCONQUE ========================== //
BSparams::BSparams(int _n, int _m, float _tmax, float _sigma, float _mu){
    n = _n;
    m = _m;
    tmax = _tmax;
    sigma = _sigma;
    mu = _mu;
}

BSparams::BSparams(){ //constructeur par défaut
    n = 0;
    m = 0;
    tmax = 0;
    sigma = 0;
    mu = 0;
}



//========================== OPTION VANILLE ==========================//
Vanilla::Vanilla(float _K, bool _call, int _n, int _m, float _tmax, float _sigma, float _mu){
    n = _n;
    m = _m;
    tmax = _tmax;
    sigma = _sigma;
    mu = _mu;
    K = _K;
    call = _call;
}

std::vector<float> Vanilla::prices(){ //prec = nombre de subdivisions
    unsigned int prec = this->n;
    std::vector<float> p;
    for(int k=0; k<int(prec*K); k++){ //entre 0 et le strike
        if(call){p.push_back(0);}
        if(not call){p.push_back( K-(float(k)/float(prec)) );}
    }
    for(int k=int(prec*K); k<int(prec); k++){
        if(call){p.push_back( float(k)/float(prec) - K);}
        if(not call){p.push_back( 0 );}
    }
    return p;
}



//========================== PORTEFEUILLE D'OPTIONS ==========================//
Portfolio::Portfolio(std::vector<Vanilla> _L, unsigned int _l){
    l = _l;
    for(unsigned int k=0; k<_l; k++){
        L.push_back(_L[k]);
    }
}
void Portfolio::add(Vanilla V){
    L.push_back(V);
    l = l+1;
}
std::vector<float> Portfolio::presentValue(unsigned int prec, int method){ //valeur du portefeuille au présent
    if(l==0){
        std::cout << "PORTEFEUILLE VIDE, VALEUR NULLE" << std::endl;
        std::vector<float> error;
        return error;
    }
    Vanilla V0 = L[0];
    unsigned int n = V0.n;
    unsigned int m = V0.m;
    std::vector<float> PresentValues(n); //matrice qui sera renvoyée

    for(unsigned int k=0; k<l; k++){
        float tmax = L[k].tmax;
        float sigma = L[k].sigma;
        float mu = L[k].mu;
        std::vector<float> prices = L[k].prices();
        BSparams par(n,m,tmax,sigma,mu);
        Matrix localSol = BSSol(par, prices, method);
        for(unsigned int v=0; v<n; v++){
            PresentValues[v] += localSol.load(v,0);
        }
    }

    for(unsigned int k=0; k<l; k++){
        PresentValues[k] /= l;
        //on divise chaque terme par l pour obtenir une valeur moyenne, qui est la valeur du portefeuille
    }

    return PresentValues;
}



//====================== ALGORITHMES DE RESOLUTION ======================//

/*La valorisation de la dérivée est une fonction du prix du sous-jacent S et
dutemps t. La fonction résolvant l'EDP de Black-Scholes doit donc prendre en
entrée les paramètres de modèle donnes ci-dessus, travaille avec une matrice
de toutes les valeurs pour tous les prix envisageables (discrets ici) et renvoie
un vecteur des valorisations pour un prix donne.*/

//création de la matrice A (voir rapport)
void matCreate(BSparams &par, Matrix &Mat){
    int m = par.m; int n = par.n;
    float dt = par.tmax/float(n);
    float a = 0.0; //variable utilisée pour remplir Mat
    float fi = 0;

    for(int i=0; i<m; i++){ //parcours des lignes
        fi = float(i);
        for(int j=0; j<n; j++){

            float eta = 1; // float eta = 1/(1-par.mu*dt);

            if( i-j<-1 or i-j>1 ){a = 0;}

            else if(i==j+1){a = 0.5*par.mu*fi*dt-0.5*par.sigma*par.sigma*fi*fi*dt;}

            else if(i==j){a = 1.0 + par.sigma*par.sigma*fi*fi*dt;}

            else if(i==j-1){a = -0.5*par.mu*fi*dt-0.5*par.sigma*par.sigma*fi*fi*dt;}

            Mat.load(i,j,eta*a);
        }
    }
}



/*YMatrix = matrice Yacine, permet de résoudre le modèle de Black Scholes sans inverser
de matrice, ce qui est plus stable empiriquement que la méthode des différences finies.*/
Matrix YMatrix(BSparams &par){
    Matrix Mat(par.n, par.m);
    float a = 0.0; //variable utilisée pour remplir Mat
    float fi = 0;
    float s = par.sigma;
    float dt = par.tmax/float(par.n);

    for(int i=0; i<par.m; i++){ //parcours des lignes
        fi = float(i);
        for(int j=0; j<par.n; j++){

            if( i-j<-1 or i-j>1 ){a = 0;}

            else if(i==j+1){a = fi*dt*(s*s-par.mu)/2;}

            else if(i==j){a = 1 - s*s*fi*dt - par.mu*dt;}

            else if(i==j-1){a = fi*dt*(s*s+par.mu)/2;}

            Mat.load(i,j,a);
        }
    }
    return Mat;
}



/* BSSol = Black-Shcholes solution : prend en entrée une instance de la classe
BSparams, définie en EqRes.hpp, ainsi que le vecteur des prix initiaux du sous-
jacent, sous forme d'un vecteur de taille par.m */
Matrix BSSol(BSparams &par, std::vector<float> &prices, unsigned int method){
    Matrix Valuations(par.m,par.n);
    Matrix Sol(par.m,1,prices);

    if(method==1){ //méthode classique des différences finies
        Matrix Mat(par.m, par.n);
        matCreate(par, Mat);
        Matrix MatInv = inverse(Mat);
        Matrix MatInv_copy = MatInv.copy(); //sinon, l'opérateur de multiplication ne fonctionne pas

        for(int t=0; t<par.n; t++){
            MatInv *= MatInv_copy;
            Matrix tempPrices = MatInv*Sol;
            for(int k=0; k<par.m; k++){
                if(tempPrices.load(k,0)<0){tempPrices.load(k,0,0);}
                Valuations.load(k,par.n-t-1,tempPrices.load(k,0));
            }
        }
    }

    else if(method ==2){//méthode sans inversion matricielle
        Matrix Mat = YMatrix(par);
        Matrix Mat_copy = Mat.copy();
        for(int t=0; t<par.n; t++){
            Mat *= Mat_copy;
            Matrix tempPrices = Mat*Sol;
            for(int k=0; k<par.m; k++){
                if(tempPrices.load(k,0)<0){tempPrices.load(k,0,0);}
                Valuations.load(k,par.n-t-1,tempPrices.load(k,0));
            }
        }
    }

    return Valuations; //retourne V à chaque instant et pour chaque S
}


Matrix BSSol(Vanilla V, unsigned int method){
    unsigned int n = V.n;
    unsigned int m = V.m;
    float tmax = V.tmax;
    float mu = V.mu;
    float sigma = V.sigma;

    std::vector<float> prices;
    prices = V.prices();

    BSparams par(n,m,tmax,sigma,mu);

    return( BSSol(par, prices, method) );
}



std::vector<float> PFvalue(Portfolio P){
    std::vector<float> res;
    for(unsigned int k=0; k<P.l; k++){
        Matrix V = BSSol(P.L[k], 2);
        res.push_back( V.load(0,k) );
    }
    for(unsigned int k=0; k<P.l; k++){
        res[k] /= float(P.l);
    }

    return res;
}



/* priceExample1 :
Générateur d'un exemple de vecteur de prix initiaux d'un sous-jacent, pour mettre
en entrée de la fonction BSSol. "gauge" correspond à la finesse du découpage de
l'intervalle [1,2], dans lequel applique une certaine fonction, en général
croissante, de la valeur du dérivé en fonction (f) de la valeur du sous-jacent.
"number" est la référence de l'exemple.

Chaque exemple correspond à une fonctionf différente :
    1. f:x -> x^2
    2. f:x -> exp(x*log(10))bhjjs
    3. f:x -> exp(1-x)
*/
std::vector<float> priceExamples(unsigned int number, unsigned int gauge){
    std::vector<float> prices(gauge);
    if(number == 1){
        for(unsigned int k=0; k<gauge; k++){
            prices[k] = 1+(float(k)/float(gauge))*(float(k)/float(gauge));
        }
    }

    else if(number == 2){
        for(unsigned int k=0; k<gauge; k++){
            prices[k] = exp(float(k)*log(10)/float(gauge-1)) - 2;
        }
    }

    else if(number == 3){
        for(unsigned int k=0; k<gauge; k++){
            prices[k] = exp(1-float(k)/float(gauge));
        }
    }
    else if(number == 4){ //call vanille
        float p = 0;
        for(unsigned int k=0; k<gauge/2; k++){prices[k] = p;}
        for(unsigned int k=gauge/2; k<gauge; k++){
            prices[k] = float(k)/gauge - 0.5;
        }
    }
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




//====================== RÉSOLUTION POUR 1 VALEUR ======================//
// (dans le cas d'un call vanille européen, voir algorithme ci-dessus
//  pour plus de détails)
float PricingExplicit(float S, float K, float r, float T, float sigma){
    float d1 = ( log(S/K) + (r + sigma*sigma /2)*T)/(sigma*sqrt(T) );
    float d2 = d1 - sigma*sqrt(T);
    float Ke = K*exp(-r*T);
    float N1 = 0.5*erfc( -d1/sqrt(2) );
    float N2 = 0.5*erfc( -d2/sqrt(2) );
    float C = S*N1 - Ke*N2;
    return C;
}


/*PricingExplicitS retourne la liste des valorisations exactes d'un call vanille européen
Elle prend en arguments le nombre de subdivisions "prec" du sous-jacent (dont la valeur
varie entre 0 et 1 par convention) et un call (de format Vanilla).
*/

std::vector<float> PricingExplicitS(unsigned int prec, Vanilla V){
    std::vector<float> C;
    for(int i=0; i<V.m ;i++){
        float uvalue = float(i)/float(prec); //=underlying value
        C.push_back( PricingExplicit(uvalue,V.K,V.mu,V.tmax,V.sigma) );
    }
    return C;
}
