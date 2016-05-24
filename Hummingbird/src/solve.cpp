/***************************************************************************
 *   Copyright (C) 2006-2011 by  Institute of Combustion Technology        *
 *   jens.henrik.goebbert@itv.rwth-aachen.de                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <iostream>
#include <math.h>

#include "data.h"
#include "metrics.h"
#include "output.h"
#include "setup.h"
#include "solve.h"

//-----------------------------------------------------
bool updateBoundaries(sData* data)
{
    // set boundary conditions for scalar s1
    double x, y, dx;
    double dhdx[2];

    for(int j = 0; j < data->nY; ++j) {
        // left boundary
        data->s1[0][j] = data->uIn * data->x[0][j];
        // right boundary
        data->s1[data->nX - 1][j] = data->uOut * data->x[data->nX - 1][j];
    }

    for(int i = 1; i < data->nX - 1; ++i) {
        // lower boundary
        // dx = (data->x[i + 1][0] - data->x[i - 1][0]) / DELD;
        dx = 1. / DELD;
        x = data->x[i][0];
        dhdx[0] = (ContourMin(x + dx) - ContourMin(x - dx)) / (2 * dx);
        data->s1[i][0] = dhdx[0] * ((data->s1[i + 1][0] - data->s1[i - 1][0]) * data->dxidx[i][0] +
                                       (4 * data->s1[i][1] - data->s1[i][2]) * data->detdx[i][0]);
        data->s1[i][0] -= (data->s1[i + 1][0] - data->s1[i - 1][0]) * data->dxidy[i][0] +
            (4 * data->s1[i][1] - data->s1[i][2]) * data->detdy[i][0];
        data->s1[i][0] /= 3 * (data->detdx[i][0] * dhdx[0] - data->detdy[i][0]);
        data->s2[i][0] = data->s2[0][0];

        // upper boundary
        // dx = (data->x[i + 1][data->nY - 1] - data->x[i - 1][data->nY - 1]) / DELD;
        x = data->x[i][data->nY - 1];
        dhdx[1] = (ContourMax(x + dx) - ContourMax(x - dx)) / (2 * dx);
        data->s1[i][data->nY - 1] =
            dhdx[1] * ((data->s1[i + 1][data->nY - 1] - data->s1[i - 1][data->nY - 1]) * data->dxidx[i][data->nY - 1] +
                          (data->s1[i][data->nY - 3] - 4 * data->s1[i][data->nY - 2]) * data->detdx[i][data->nY - 1]);
        data->s1[i][data->nY - 1] -=
            (data->s1[i + 1][data->nY - 1] - data->s1[i - 1][data->nY - 1]) * data->dxidy[i][data->nY - 1] +
            (data->s1[i][data->nY - 3] - 4 * data->s1[i][data->nY - 2]) * data->detdy[i][data->nY - 1];
        data->s1[i][data->nY - 1] /= 3 * (data->detdy[i][data->nY - 1] - data->detdx[i][data->nY - 1] * dhdx[1]);
        data->s2[i][data->nY - 1] = data->s2[0][data->nY - 1];
    }

    // printValue(data, "s boundaries", data->s1);
    return true;
}

bool solve(sData* data)
{
    std::cout << "\nSolve:\n-------\n";

    if(data->solverType == GAUSSSEIDEL) {
        if(!(gaussseidel(data, data->s1) && gaussseidel(data, data->s2))) {
            return false;
        }
    } else if(data->solverType == JACOBI) {
        if(!(jacobi(data, data->s1) && jacobi(data, data->s2))) {
            return false;
        }
    } else if(data->solverType == THOMAS) {
        if(!(thomas(data, data->s1) && thomas(data, data->s2))) {
            return false;
        }
    }

    return true;
}

//--- Gauss-Seidel solver -----------------------------
bool gaussseidel(sData* data, double** s)
{
    int curIter = 0;
    double curResidual = MAXDOUBLE;
    double curMaxResidual = MAXDOUBLE;
    double sBar;

    std::cout << "\n\r\tRunning Gauss-Seidel solver... ";

    while(curIter < data->maxIter && ABS(curMaxResidual) > data->maxResidual) {
        curIter++;

        curMaxResidual = 0;

        for(int i = 1; i < data->nX - 1; i++) {
            for(int j = 1; j < data->nY - 1; j++) {

                curResidual = -s[i][j];

                sBar = (s[i + 1][j] + s[i - 1][j]) * data->a1[i][j] + (s[i + 1][j] - s[i - 1][j]) * data->a4[i][j] / 2 +
                    (s[i][j + 1] + s[i][j - 1]) * data->a2[i][j] + (s[i][j + 1] - s[i][j - 1]) * data->a5[i][j] / 2;
                sBar += (s[i + 1][j + 1] - s[i - 1][j + 1] + s[i - 1][j - 1] - s[i + 1][j - 1]) * data->a3[i][j] / 4;
                sBar /= 2 * (data->a1[i][j] + data->a2[i][j]);

                s[i][j] = s[i][j] + data->facRelax * (sBar - s[i][j]);

                curResidual += s[i][j];

                if(ABS(curResidual) > ABS(curMaxResidual))
                    curMaxResidual = ABS(curResidual);
            }
        }

        if(data->potentialFunc == CHANNELFLOW) {
            updateBoundaries(data);
        }
    }

    std::cout << "done." << std::endl;
    std::cout << "\n\tGauss-Seidel iterations: " << curIter << std::endl;
    std::cout << "\n\tMaximum residuum       : " << curMaxResidual << std::endl;

    return true;
}

//--- Jacobi solver -----------------------------------
bool jacobi(sData* data, double** s)
{
    int curIter = 0;
    double curResidual = MAXDOUBLE;
    double curMaxResidual = MAXDOUBLE;
    double sBar;
    double** sOld = new double*[data->nX];

    for(int i = 0; i < data->nX; ++i) {
        sOld[i] = new double[data->nY];

        for(int j = 0; j < data->nY; ++j)
            sOld[i][j] = s[i][j];
    }

    std::cout << "\r\tRunning Jacobi solver... ";

    while(curIter < data->maxIter && ABS(curMaxResidual) > data->maxResidual) {
        curIter++;

        curMaxResidual = 0;

        for(int i = 1; i < data->nX - 1; i++) {
            for(int j = 1; j < data->nY - 1; j++) {

                curResidual = -sOld[i][j];

                sBar = (sOld[i + 1][j] + sOld[i - 1][j]) * data->a1[i][j] +
                    (sOld[i + 1][j] - sOld[i - 1][j]) * data->a4[i][j] / 2 +
                    (sOld[i][j + 1] + sOld[i][j - 1]) * data->a2[i][j] +
                    (sOld[i][j + 1] - sOld[i][j - 1]) * data->a5[i][j] / 2;
                sBar += (sOld[i + 1][j + 1] - sOld[i - 1][j + 1] + sOld[i - 1][j - 1] - sOld[i + 1][j - 1]) *
                    data->a3[i][j] / 4;
                sBar /= 2 * (data->a1[i][j] + data->a2[i][j]);

                s[i][j] = sOld[i][j] + data->facRelax * (sBar - sOld[i][j]);

                curResidual += s[i][j];

                if(ABS(curResidual) > ABS(curMaxResidual))
                    curMaxResidual = ABS(curResidual);
            }
        }

        if(data->potentialFunc == CHANNELFLOW) {
            updateBoundaries(data);
        }

        for(int i = 0; i < data->nX; ++i) {
            for(int j = 0; j < data->nY; ++j)
                sOld[i][j] = s[i][j];
        }
    }

    std::cout << "done." << std::endl;
    std::cout << "\n\tJacobi iterations: " << curIter << std::endl;
    std::cout << "\n\tMaximum residuum       : " << curMaxResidual << std::endl;

    return true;
}

//--- Thomas algorithm --------------------------------
bool thomas(sData* data, double** s)
{
    int curIter = 0;
    double curResidual = MAXDOUBLE;
    double curMaxResidual = MAXDOUBLE;
    double sBar;
    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    double* cPrimeX = new double[data->nX - 2];
    double* dPrimeX = new double[data->nX - 2];
    double* cPrimeY = new double[data->nY - 2];
    double* dPrimeY = new double[data->nY - 2];

    std::cout << "\r\tRunning Thomas solver... ";

    while(curIter < data->maxIter && ABS(curMaxResidual) > data->maxResidual) {
        curIter++;

        curMaxResidual = 0;

        if(true || curIter % 2 == 0) {

            for(int j = 1; j < data->nY - 1; ++j) {
                b = 2 * (data->a1[1][j] + data->a2[1][j]);
                c = -(data->a1[1][j] + data->a4[1][j] / 2);
                d = s[0][j] * (data->a1[1][j] - data->a4[1][j] / 2) + data->a2[1][j] * (s[1][j + 1] + s[1][j - 1]) +
                    data->a3[1][j] / 4 * (s[2][j + 1] - s[0][j + 1] - s[2][j - 1] + s[0][j - 1]) +
                    data->a5[1][j] / 2 * (s[1][j + 1] - s[1][j - 1]);
                cPrimeX[0] = c / b;
                dPrimeX[0] = d / b;
                for(int i = 2; i < data->nX - 2; ++i) {
                    a = data->a4[i][j] / 2 - data->a1[i][j];
                    b = 2 * (data->a1[i][j] + data->a2[i][j]);
                    c = -(data->a1[i][j] + data->a4[i][j] / 2);
                    d = data->a2[i][j] * (s[i][j + 1] + s[i][j - 1]) +
                        data->a3[i][j] / 4 * (s[i + 1][j + 1] - s[i - 1][j + 1] - s[i + 1][j - 1] + s[i - 1][j - 1]) +
                        data->a5[i][j] / 2 * (s[i][j + 1] - s[i][j - 1]);
                    cPrimeX[i - 1] = c / (b - cPrimeX[i - 2] * a);
                    dPrimeX[i - 1] = (d - dPrimeX[i - 2] * a) / (b - cPrimeX[i - 2] * a);
                }
                a = data->a4[data->nX - 2][j] / 2 - data->a1[data->nX - 2][j];
                b = 2 * (data->a1[data->nX - 2][j] + data->a2[data->nX - 2][j]);
                d = s[data->nX - 1][j] * (data->a1[data->nX - 2][j] + data->a4[data->nX - 2][j] / 2) +
                    data->a2[data->nX - 2][j] * (s[data->nX - 2][j + 1] + s[data->nX - 2][j - 1]) +
                    data->a3[data->nX - 2][j] / 4 * (s[data->nX - 1][j + 1] - s[data->nX - 3][j + 1] -
                                                        s[data->nX - 1][j - 1] + s[data->nX - 3][j - 1]) +
                    data->a5[data->nX - 2][j] / 2 * (s[data->nX - 2][j + 1] - s[data->nX - 2][j - 1]);
                cPrimeX[data->nX - 3] = c / (b - cPrimeX[data->nX - 4] * a);
                dPrimeX[data->nX - 3] = (d - dPrimeX[data->nX - 4] * a) / (b - cPrimeX[data->nX - 4] * a);
                curResidual = -s[data->nX - 2][j];
                sBar = dPrimeX[data->nX - 3];
                s[data->nX - 2][j] = s[data->nX - 2][j] + data->facRelax * (sBar - s[data->nX - 2][j]);
                curResidual += s[data->nX - 2][j];
                if(ABS(curResidual) > ABS(curMaxResidual))
                    curMaxResidual = ABS(curResidual);
                for(int i = data->nX - 3; i > 0; --i) {
                    curResidual = -s[i][j];
                    sBar = dPrimeX[i - 1] - cPrimeX[i - 1] * s[i + 1][j];
                    s[i][j] = s[i][j] + data->facRelax * (sBar - s[i][j]);
                    curResidual += s[i][j];
                    if(ABS(curResidual) > ABS(curMaxResidual))
                        curMaxResidual = ABS(curResidual);
                }
            }
        } else {
            for(int i = 1; i < data->nX - 1; ++i) {
                b = 2 * (data->a1[i][1] + data->a2[i][1]);
                c = -(data->a2[i][1] + data->a5[i][1] / 2);
                d = s[i][0] * (data->a2[i][1] - data->a5[i][1] / 2) + data->a1[i][1] * (s[i + 1][1] + s[i - 1][1]) +
                    data->a3[i][1] / 4 * (s[i + 1][2] - s[i - 1][2] - s[i + 1][0] + s[i - 1][0]) +
                    data->a4[i][1] / 2 * (s[i + 1][1] - s[i - 1][1]);
                cPrimeY[0] = c / b;
                dPrimeY[0] = d / b;
                for(int j = 2; j < data->nY - 2; ++j) {
                    a = data->a5[i][j] / 2 - data->a2[i][j];
                    b = 2 * (data->a1[i][j] + data->a2[i][j]);
                    c = -(data->a2[i][j] + data->a5[i][j] / 2);
                    d = data->a1[i][j] * (s[i + 1][j] + s[i - 1][j]) +
                        data->a3[i][j] / 4 * (s[i + 1][j + 1] - s[i - 1][j + 1] - s[i + 1][j - 1] + s[i - 1][j - 1]) +
                        data->a4[i][j] / 2 * (s[i + 1][j] - s[i - 1][j]);
                    cPrimeY[j - 1] = c / (b - cPrimeY[j - 2] * a);
                    dPrimeY[j - 1] = (d - dPrimeY[j - 2] * a) / (b - cPrimeY[j - 2] * a);
                }
                a = data->a5[i][data->nY - 2] / 2 - data->a2[i][data->nY - 2];
                b = 2 * (data->a1[i][data->nY - 2] + data->a2[i][data->nY - 2]);
                d = s[i][data->nY - 1] * (data->a2[i][data->nY - 2] + data->a5[i][data->nY - 2] / 2) +
                    data->a1[i][data->nY - 2] * (s[i + 1][data->nY - 2] + s[i - 1][data->nY - 2]) +
                    data->a3[i][data->nY - 2] / 4 * (s[i + 1][data->nY - 1] - s[i - 1][data->nY - 1] -
                                                        s[i + 1][data->nY - 3] + s[i - 1][data->nY - 3]) +
                    data->a4[i][data->nY - 2] / 2 * (s[i + 1][data->nY - 2] - s[i - 1][data->nY - 2]);
                cPrimeY[data->nY - 3] = c / (b - cPrimeY[data->nY - 4] * a);
                dPrimeY[data->nY - 3] = (d - dPrimeY[data->nY - 4] * a) / (b - cPrimeY[data->nY - 4] * a);
                curResidual = -s[i][data->nY - 2];
                sBar = dPrimeY[data->nY - 3];
                s[i][data->nY - 2] = s[i][data->nY - 2] + data->facRelax * (sBar - s[i][data->nY - 2]);
                curResidual += s[i][data->nY - 2];
                if(ABS(curResidual) > ABS(curMaxResidual))
                    curMaxResidual = ABS(curResidual);
                for(int j = data->nY - 3; j > 0; --j) {
                    curResidual = -s[i][j];
                    sBar = dPrimeY[j - 1] - cPrimeY[j - 1] * s[i][j + 1];
                    s[i][j] = s[i][j] + data->facRelax * (sBar - s[i][j]);
                    curResidual += s[i][j];
                    if(ABS(curResidual) > ABS(curMaxResidual))
                        curMaxResidual = ABS(curResidual);
                }
            }
        }
        if(data->potentialFunc == CHANNELFLOW) {
            updateBoundaries(data);
        }
    }

    std::cout << "done." << std::endl;
    std::cout << "\n\tThomas iterations: " << curIter << std::endl;
    std::cout << "\n\tMaximum residuum       : " << curMaxResidual << std::endl;

    return true;
}
