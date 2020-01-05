#include <QtCharts/QLineSeries>
#include "eqres.h"

#ifndef DISPLAY_H
#define DISPLAY_H

using namespace QtCharts;

class BSparams;

QLineSeries * ResExample(bool time, int t, unsigned int prec,
                         BSparams par, unsigned int exampleNumber,unsigned int method);

int displayTest(int argc, char *argv[],bool time, int tranche,
                unsigned int prec, BSparams par, unsigned int exampleNumber, unsigned int method);

int displayVector(std::vector<float> X, std::vector<std::vector<float>> Y, int argc, char *argv[]);

#endif // DISPLAY_H
