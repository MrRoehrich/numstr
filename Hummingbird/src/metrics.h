#ifndef METRICS_H
#define METRICS_H
#include "data.h"


// functions
double Xi(sData* data, double x, double y);
double Et(sData* data, double x, double y);
double X(sData* data, double xi, double et);
double Y(sData* data, double xi, double et);
double ContourMax(double x);
double ContourMin(double x);

#endif
