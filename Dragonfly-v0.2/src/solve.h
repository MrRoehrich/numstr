/***************************************************************************
 *   Copyright (C) 2006-2011 by  Institute of Combustion Technology        *
 *   d.mayer@itv.rwth-aachen.de                                            *
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
#ifndef CALC_H
#define CALC_H

#include "data.h"

bool solveCalcFlux(sData* data);
bool solvePe(sData* data);
void calcFluxCentral(sData* data);
void calcFluxUpwind(sData* data);
double A(double Pe);
double calcPe(sData* data, double velocity, double length);
void calcVelocityField(sData* data, double deltaT);
void calcScalarField(sData* data, double deltaT);
bool solveSimple(sData* data);
bool solveSimpler(sData* data);
void calcCoeff(sData* data,
    double diffCoef,
    double deltaT,
    double dx,
    double dy,
    double vn,
    double ue,
    double vs,
    double uw,
    double& an,
    double& ae,
    double& as,
    double& aw,
    double& ap,
    double& apTilde);

#endif
