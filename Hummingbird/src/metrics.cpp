#include <iostream>
#include <math.h>

#include "data.h"
#include "metrics.h"

//---- countur functions ------------------------------
double ContourMax(double x)
{
    double cmax;
     cmax = 100.;
    // cmax = 1. +(3.-1.)/(1.-0.) * (x-0.);
    // cmax = atanl(x) + 1.;
    // cmax = atanl(x) + 2.;
    // cmax = atanl(x) + 8.;
    return cmax;
}

double ContourMin(double x)
{
    double cmin;
    cmin = 0.01;
    // cmin = -1 +(-3.-(-1))/(30.-0) * (x-0);
    // cmin = 0.;
    // cmin = atanl(x) - 2.;
    return cmin;
}

//---- coordinate transformation ----------------------
double Xi(sData* data, double x, double y)
{
    double xi = (data->nX - 1) / (data->xMax - data->xMin) * (x - data->xMin);
    return xi;
}

double Et(sData* data, double x, double y)
{
    double et = (data->nY - 1) / (ContourMax(x) - ContourMin(x)) * (y - ContourMin(x));
    return et;
}

double X(sData* data, double xi, double et)
{
    double x = data->xMin + (data->xMax - data->xMin) * xi / (data->nX - 1);
    return x;
}

double Y(sData* data, double xi, double et)
{
    double x = X(data, xi, et);
    double y = ContourMin(x) + (ContourMax(x) - ContourMin(x)) * et / (data->nY - 1);
    return y;
}
