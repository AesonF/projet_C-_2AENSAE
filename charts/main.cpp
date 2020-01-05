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



/*std::vector<float> portCreator(){
    std::vector<Vanilla> Van;
    Vanilla V1(0.6,true,50,50,1,0.1,0.1);
    Vanilla V2(0.4,false,50,50,1,0.1,0.1);
    Van.push_back(V1);
    Van.push_back(V2);
    Portfolio P(Van,2);
    std::vector<float> Values = P.presentValue(50,2);
    return Values;
}*/


/*void genericPlotTest(int argc, char *argv[]){
    std::vector<float> X,Y1,Y2;
    for(int k=0; k<100; k++){
        X.push_back(float(k));
        Y1.push_back(float(k*k));
        Y2.push_back(float(k*k/2));
    }
    std::vector<std::vector<float>> Y;
    Y.push_back(Y1);
    Y.push_back(Y2);s
    displayVector(X,Y,argc, argv);
}*/



int main(int argc, char *argv[]){
    unsigned int prec = 100;
    BSparams par(prec, prec,0.5,0.01,0.3);
    int b = displayTest(argc, argv, true, 0, prec, par, 2, 2);
    return b;

    /*std::vector<float> X;
    for(int k=0; k<50; k++){
        X.push_back(float(k)/50);
    }
    std::vector<std::vector<float>> Y;
    Y.push_back(portCreator());
    displayVector(X, Y,argc,argv);*/
}
