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

class Matrix;



/* Resexample:
 * Cette fonction facilite l'affichage avec les fonctions du package QtGraphs de Qt.
 * En effet, on affiche un graphe avec des listes d'abscisses et d'ordonnées sous la
 * forme d'une variable de type "series". ResExample prend en entrée trois arguments :
 *      - "time" : true si on regarde les valeurs V(.,t) pour un temps t donné, faux si on
 *        regarde les valeurs V(s,.) pour une valeur s du sous-jacent donnée.
 *      - "tranche" : spécifie la valeur fixe de t (si time=true) ou de s (sinon)
 *      - "prec" : nombre de subdivisions de s et de t dans l'approximation de V(s,t)
 *        ("prec" correspond à Delta_t dans le rapport)
 */
QLineSeries * ResExample(bool time, int tranche, unsigned int prec,
                         BSparams par, unsigned int exampleNumber, unsigned int method) {

    std::vector<float> prices = priceExamples(exampleNumber, prec);
    Matrix V = BSSol(par, prices, method); //toutes les valorisations

    QLineSeries *series = new QLineSeries;

    if(time==true){ //on fixe le temps
        for(unsigned int k=0; k<prec; k++){
            float r = V.load(k,tranche);
            qreal yval = qreal(r);
            series->append(float(k)/float(prec),yval); //S varie entre 0 et 1
        }
    }
    else{ //on fixe la valeur du dérivé à maturité
        for(unsigned int k=0; k<prec; k++){
            float r = V.load(tranche, k);
            qreal yval = qreal(r);
            series->append(float(k*par.tmax)/float(prec),yval); //S varie entre 0 et 1
        }
    }
    return series;
}



int displayTest(int argc, char *argv[],bool time, int tranche,
                unsigned int prec, BSparams par, unsigned int exampleNumber, unsigned int method) {

    // 1. liste des valorisations selon l'algorithme de résolution (méthodes 1 et 2)
    QLineSeries * series = ResExample(time, tranche , prec, par, exampleNumber, method);

    // 2. liste des valorisations à maturité
    QLineSeries *Pseries = new QLineSeries;
    std::vector<float> prices = priceExamples(exampleNumber, prec);
    for(unsigned int k=0; k<prec; k++){
        qreal yval = qreal(prices[k]);
        Pseries->append(float(k)/prec,yval);
    }

    QApplication a(argc, argv);

    QChart *chart = new QChart;
    chart->legend()->show();

    chart->addSeries(series);
    if(time==true){chart->addSeries(Pseries);}

    // 3. liste des valorisations exactes d'un call vanille européen
    if(time==true and exampleNumber==4){
        float K = 0.5;
        bool call = true;
        Vanilla V = Vanilla(K,call,par.n,par.m,par.tmax,par.sigma,par.mu);
        std::vector<float> evalues = PricingExplicitS(prec, V);
        QLineSeries *Tseries = new QLineSeries;
        for(unsigned int k=0; k<prec; k++){
            qreal yval = qreal(evalues[k]);
            Tseries->append(float(k)/prec,yval);
        }
        chart->addSeries(Tseries);
    }

    chart->createDefaultAxes();

    QFont font;
    font.setPixelSize(18);
    chart->setTitleFont(font);
    if(time==true){
        chart->setTitle("Value of derivative at two instants, as a function of underlying value");
    }
    else{
        chart->setTitle("Value of derivative over time, with fixed price at maturity");
    }
    QPen pen(QRgb(0x0000000));
    pen.setWidth(2);
    series->setPen(pen);

    //chart->setAnimationOptions(QChart::AllAnimations);

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(chartView);
    window.resize(700,500);
    window.show();

    return a.exec();
}



int displayVector(std::vector<float> X, std::vector<std::vector<float>> Y, int argc, char *argv[]){
    QApplication a(argc, argv);

    QChart *chart = new QChart;
    chart->legend()->show();

    QPen pen(QRgb(0x0000000));
    pen.setWidth(2);

    for(unsigned int r=0; r<Y.size(); r++){
        QLineSeries * series = new QLineSeries;
        for(unsigned int k=0; k<X.size(); k++){
            float y = Y[r][k];
            qreal xval = qreal(X[k]);
            qreal yval = qreal(y);
            series->append(xval,yval);
        }
        chart->addSeries(series);
        series->setPen(pen);
    }

    chart->createDefaultAxes();

    QFont font;
    font.setPixelSize(18);
    chart->setTitleFont(font);
    chart->setTitle("Generic plot");

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(chartView);
    window.resize(700,500);
    window.show();

    return a.exec();
}
