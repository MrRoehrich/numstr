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
#ifndef DATA_H
#define DATA_H

#include <limits>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ABS(x) (((x) > 0) ? (x) : -(x))

#define MAXDOUBLE (std::numeric_limits<double>::max())
#define MINDOUBLE (std::numeric_limits<double>::min())
#define MAXFLOAT (std::numeric_limits<float>::max())
#define MINFLOAT (std::numeric_limits<float>::min())
#define MAXINT (std::numeric_limits<int>::max())
#define MININT (std::numeric_limits<int>::min())

#define DELD 1.e6
#define DELDD 1.e3

#define PI 3.14

#define CHANNELFLOW -1
#define PARALLELFLOW 0
#define STAGNATION_POINT 1
#define SOURCE 2
#define POTENTIAL_VORTEX 3
#define DIPOL 4

#define GAUSSSEIDEL 0
#define JACOBI 1
#define THOMAS 2

struct sData {
    sData() {};

    // grid & geometry settings
    int nX;
    int nY;
    double xMin;
    double xMax;

    // potential flow settings
    int potentialFunc;
    double uIn;
    double uOut;
    double uInfty;
    double vInfty;
    double magnitude;

    // solver settings
    int solverType;
    int maxIter;
    double maxResidual;
    double facRelax;

    // fields
    double** x;
    double** y;
    double** xi;
    double** et;
    double** dxidx;
    double** dxidy;
    double** detdx;
    double** detdy;
    double** ddxidx;
    double** ddxidy;
    double** ddetdx;
    double** ddetdy;
    double** a1;
    double** a2;
    double** a3;
    double** a4;
    double** a5;
    double** s1;  // field for potentia lfunction
    double** s2;  // field for stream function
    double** ds1dxi;
    double** ds1det;
    double** u;
    double** v;
    double** cp;
};

double** allocGrid1Mem(const sData* const data, const double preset);
void freeGrid1Mem(const sData* const data, double** mem);
double*** allocGridXMem(const sData* const data, const int vSize, const double preset);

#endif
