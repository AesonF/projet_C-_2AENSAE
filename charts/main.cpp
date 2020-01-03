
//Standard packages:
#include <string>
#include <vector>
#include <cmath>

//Qt packages:
#include <QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QCategoryAxis>

//Qt headers:
#include "mainwindow.h"

//My headers:
#include "matrix.h"
#include "eqres.h"
#include "display.h"

QT_CHARTS_USE_NAMESPACE

/* EXEMPLE DE CRÉATION D'UN PORTEFEUILLE
void PortCreator(){
    std::vector<Vanilla> Van;
    Vanilla V1(0.6,true,50,50,1,0.1,0.1);
    Vanilla V2(0.4,false,50,50,1,0.1,0.1);
    Van.push_back(V1);
    Van.push_back(V2);
    Portfolio P(Van,2);
    Matrix Values = P.presentValue(50,2);
    Values.show();
}
*/

int main(int argc, char *argv[]){
    unsigned int prec = 100;
    BSparams par(prec, prec,0.5,0.01,0.3);
    int b = displayTest(argc, argv, true, 0, prec, par, 4, 1);
    return b;

    /*float n = PricingExplicit(3850,4100,0.0125,1,0.0168);
    std::cout << "Explicit price: " << n << std::endl;*/
}
