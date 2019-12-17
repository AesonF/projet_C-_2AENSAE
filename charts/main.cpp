//Standard packages:
#include <string>
#include <vector>

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

QT_CHARTS_USE_NAMESPACE

class Matrix;

std::vector<float> ResExample() {
    unsigned int prec = 100;
    std::vector<float> prices = priceExample2(prec);
    BSparams par(prec,prec,10.0,0.1,0.02);
    Matrix V = BSSol(par, prices);
    std::vector<float> res;
    for(int k=0; k<prec; k++){
        res.push_back(V.load(100,k));
    }

    return res;
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    std::vector<float> prices_at_maturity = ResExample();

    QLineSeries *series = new QLineSeries;
    for(int k=0; k<100; k++){
        qreal yval = qreal(prices_at_maturity[k]);
        series->append(k,yval);
    }

    QChart *chart = new QChart;
    chart->legend()->show();
    chart->addSeries(series);
    chart->createDefaultAxes();

    QFont font;
    font.setPixelSize(18);
    chart->setTitleFont(font);
    chart->setTitle("Approximated values over time");
    QPen pen(QRgb(0x0000000));
    pen.setWidth(2);
    series->setPen(pen);

    //chart->setAnimationOptions(QChart::AllAnimations);

    QCategoryAxis *axisX = new QCategoryAxis();
    for(int k=0; k<100; k++){
        axisX->append("",k*k);
    }

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(chartView);
    window.resize(800,800);
    window.show();

    return a.exec();
}
