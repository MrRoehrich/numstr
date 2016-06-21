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
#include <iostream>

#include <math.h>
#include <stdio.h>

#include "data.h"
#include "output.h"
#include "solve.h"

//------------------------------------------------------
double A(double Pe, sData* data)
{
    if(data->solverType == CENTRAL)
        return 1. - ABS(Pe) / 2.;
    else if(data->solverType == UPWIND)
        return 1.;
    else if(data->solverType == HYBRID)
        return MAX(0., 1. - ABS(Pe) / 2.);
    else if(data->solverType == POWER)
        return MAX(0., powf(1. - ABS(Pe) / 10., 5.));
    else if(data->solverType == EXPONENTIAL)
        if(data->alpha == 0)
            return 1.;
    return ABS(Pe) / (exp(ABS(Pe)) - 1.);
}

//------------------------------------------------------
double calcPe(sData* data, double velocity, double length)
{
    if(data->alpha > EPS) {
        return data->rho * velocity / data->alpha * length;
    } else {
        if(velocity > 0.)
            return MAXDOUBLE;
        else
            return -MAXDOUBLE;
    }
}

//------------------------------------------------------
bool solveCalcFlux(sData* data)
{
    std::cout << "\nCalculation:\n------------\n";
    static sCell* curCell = 0;
    static sFace* ym = 0;
    static sFace* yp = 0;
    static sFace* xm = 0;
    static sFace* xp = 0;

    double curTime = 0.;
    int curIter = 0.;
    double deltaT = data->maxTime / data->maxIter;

    // write initial conditions
    std::cout << "\nOutput... " << curIter << "\n";
    if(!output(data, curIter)) {
        std::cout << "ERROR while data output...exiting";
        getchar();
        return 1;
    }

    while(curTime < data->maxTime && curIter < data->maxIter) {

        curIter++;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        calcVelocityField(data, deltaT);

        if(data->solverType == CENTRAL)
            calcFluxCentral(data);
        else if(data->solverType == UPWIND)
            calcFluxUpwind(data);

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if(curCell->bTypeScalar == 1)
                continue;
            ym = curCell->faces[YM];
            yp = curCell->faces[YP];
            xm = curCell->faces[XM];
            xp = curCell->faces[XP];

            // fluxes f* (horizontally) and g* (vertically)
            curCell->fluxBalance = ym->numFlux[M] * ym->dy - ym->numFlux[P] * ym->dx;
            curCell->fluxBalance += -yp->numFlux[M] * yp->dy + yp->numFlux[P] * yp->dx;
            curCell->fluxBalance += -xm->numFlux[M] * xm->dy - xm->numFlux[P] * xm->dx;
            curCell->fluxBalance += xp->numFlux[M] * xp->dy + xp->numFlux[P] * xp->dx;

            curCell->phi += deltaT / data->rho * curCell->s;
            curCell->phi -= deltaT / (curCell->volume * data->rho) * curCell->fluxBalance;
        }

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if((cId + data->nCellsX) % data->nCellsX == 0) { // linker Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[XM]->u = curCell->faces[XP]->u;
                    curCell->faces[XM]->v = curCell->faces[XP]->v;
                }
            } else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) { // rechter Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[XP]->u = curCell->faces[XM]->u;
                    curCell->faces[XP]->v = curCell->faces[XM]->v;
                }
            } else if(cId < data->nCellsX) { // unterer Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[YM]->u = curCell->faces[YP]->u;
                    curCell->faces[YM]->v = curCell->faces[YP]->v;
                }
            } else if(data->nCells - (cId + 1) < data->nCellsX) { // oberer Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[YP]->u = curCell->faces[YM]->u;
                    curCell->faces[YP]->v = curCell->faces[YM]->v;
                }
            }
        }

        curTime += deltaT;
        // write output
        std::cout << "Output... " << curIter << "\n";
        if(!output(data, curIter)) {
            std::cout << "ERROR while data output...exiting";
            getchar();
            return 1;
        }
    }
    return true;
}

bool solvePe(sData* data)
{
    std::cout << "\nCalculation:\n------------\n";
    static sCell* curCell = 0;

    double apTilde = 0.;
    double ap = 0.;
    double an = 0.;
    double ae = 0.;
    double as = 0.;
    double aw = 0.;
    double b = 0.;

    double D = 0.;
    double Pe = 0.;
    double f = 0.;
    double g = 0.;

    double curTime = 0.;
    int curIter = 0.;
    double deltaT = data->maxTime / data->maxIter;

    // write initial conditions
    std::cout << "Output... " << curIter << "\n";
    if(!output(data, curIter)) {
        std::cout << "ERROR while data output...exiting";
        getchar();
        return 1;
    }

    while(curTime < data->maxTime && curIter < data->maxIter) {
        curIter++;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        calcVelocityField(data, deltaT);

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if(curCell->bTypeScalar != 0)
                continue;

            // an
            D = data->alpha / curCell->faces[YP]->dx;
            Pe = calcPe(data, curCell->faces[YP]->v, curCell->faces[YP]->dx);
            g = data->rho * curCell->faces[YP]->v;
            an = D * curCell->faces[YP]->dx * A(Pe, data) + MAX(-g * curCell->faces[YP]->dx, 0.);
            // ae
            D = data->alpha / curCell->faces[XP]->dy;
            Pe = calcPe(data, curCell->faces[XP]->u, curCell->faces[XP]->dy);
            f = data->rho * curCell->faces[XP]->u;
            ae = D * curCell->faces[XP]->dy * A(Pe, data) + MAX(-f * curCell->faces[XP]->dy, 0.);
            // as
            D = data->alpha / curCell->faces[YM]->dx;
            Pe = calcPe(data, curCell->faces[YM]->v, curCell->faces[YM]->dx);
            g = data->rho * curCell->faces[YM]->v;
            as = D * curCell->faces[YM]->dx * A(Pe, data) + MAX(g * curCell->faces[YM]->dx, 0.);
            // aw
            D = data->alpha / curCell->faces[XM]->dy;
            Pe = calcPe(data, curCell->faces[XM]->u, curCell->faces[XM]->dy);
            f = data->rho * curCell->faces[XM]->u;
            aw = D * curCell->faces[XM]->dy * A(Pe, data) + MAX(f * curCell->faces[XM]->dy, 0.);
            // ap
            ap = data->rho * curCell->faces[XM]->dy * curCell->faces[YM]->dx / deltaT;
            // b
            b = ap * curCell->phi;
            // apTilde
            apTilde = an + ae + as + aw + ap;

            curCell->phi = an * curCell->neighCells[YP]->phi + ae * curCell->neighCells[XP]->phi +
                as * curCell->neighCells[YM]->phi + aw * curCell->neighCells[XM]->phi + b;
            curCell->phi /= apTilde;
        }

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if((cId + data->nCellsX) % data->nCellsX == 0) { // linker Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[XM]->u = curCell->faces[XP]->u;
                    curCell->faces[XM]->v = curCell->faces[XP]->v;
                }
            } else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) { // rechter Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[XP]->u = curCell->faces[XM]->u;
                    curCell->faces[XP]->v = curCell->faces[XM]->v;
                }
            } else if(cId < data->nCellsX) { // unterer Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[YM]->u = curCell->faces[YP]->u;
                    curCell->faces[YM]->v = curCell->faces[YP]->v;
                }
            } else if(data->nCells - (cId + 1) < data->nCellsX) { // oberer Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[YP]->u = curCell->faces[YM]->u;
                    curCell->faces[YP]->v = curCell->faces[YM]->v;
                }
            }
        }

        curTime += deltaT;
        // write output
        std::cout << "Output... " << curIter << "\n";
        if(!output(data, curIter)) {
            std::cout << "ERROR while data output...exiting";
            getchar();
            return 1;
        }
    }

    std::cout << "\n";
    return true;
}

//------------------------------------------------------
void calcFluxCentral(sData* data)
{
    static sFace* curFace = 0;

    double conv, diff, dx, dy;
    // compute numerical flux of each face
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(curFace->bType == 0) {

            // convective part of f*
            conv = data->rho * curFace->u * (curFace->neighCells[P]->phi + curFace->neighCells[M]->phi) / 2.;
            // diffusive part of f*
            dx = curFace->neighCells[P]->x - curFace->neighCells[M]->x;
            if(dx != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dx;
            // f* = f*(conv) + f*(diff)
            curFace->numFlux[M] = conv + diff;

            // convective part of g*
            conv = data->rho * curFace->v * (curFace->neighCells[P]->phi + curFace->neighCells[M]->phi) / 2.;
            // diffusive part of g*
            dy = curFace->neighCells[P]->y - curFace->neighCells[M]->y;
            if(dy != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dy;
            // g* = g*(conv) + g*(diff)
            curFace->numFlux[P] = conv + diff;
        } else if(curFace->bType == 2) {
            curFace->numFlux[M] = curFace->bValueX;
            curFace->numFlux[P] = curFace->bValueY;
        }
    }
}

//------------------------------------------------------
void calcFluxUpwind(sData* data)
{
    static sFace* curFace = 0;

    double conv, diff, dx, dy;
    // compute numerical flux of each face
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(curFace->bType == 0) {

            // convective part of f*
            if(curFace->u >= 0)
                conv = data->rho * curFace->u * curFace->neighCells[M]->phi;
            else
                conv = data->rho * curFace->u * curFace->neighCells[P]->phi;
            // diffusive part of f*
            dx = curFace->neighCells[P]->x - curFace->neighCells[M]->x;
            if(dx != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dx;
            // f* = f*(conv) + f*(diff)
            curFace->numFlux[M] = conv + diff;

            // convective part of g*
            if(curFace->v >= 0)
                conv = data->rho * curFace->v * curFace->neighCells[M]->phi;
            else
                conv = data->rho * curFace->v * curFace->neighCells[P]->phi;
            // diffusive part of g*
            dy = curFace->neighCells[P]->y - curFace->neighCells[M]->y;
            if(dy != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dy;
            // g* = g*(conv) + g*(diff)
            curFace->numFlux[P] = conv + diff;
        } else if(curFace->bType == 2) {
            curFace->numFlux[M] = curFace->bValueX;
            curFace->numFlux[P] = curFace->bValueY;
        }
    }
}

//------------------------------------------------------
void calcVelocityField(sData* data, double deltaT)
{
    sCell* curCell = 0;
    int cId;

    double apTilde = 0.;
    double ap = 0.;
    double an = 0.;
    double ae = 0.;
    double as = 0.;
    double aw = 0.;
    double b = 0.;

    double D = 0.;
    double Pe = 0.;
    double f = 0.;
    double g = 0.;

    for(int nY = 2; nY < data->nCellsY; ++nY) {
        for(int nX = 2; nX < data->nCellsX; ++nX) {
            cId = (nX - 1) + (nY - 1) * data->nCellsX;
            curCell = &data->cells[cId];

            // u
            // an
            D = data->alpha / (curCell->neighCells[XP]->x - curCell->x);
            Pe = calcPe(data, curCell->neighCells[YP]->faces[XP]->v, (curCell->neighCells[XP]->x - curCell->x));
            g = data->rho * curCell->neighCells[YP]->faces[XP]->v;
            an = D * (curCell->neighCells[XP]->x - curCell->x) * A(Pe, data) +
                MAX(-g * (curCell->neighCells[XP]->x - curCell->x), 0.);
            // ae
            D = data->alpha / curCell->faces[XP]->dy;
            Pe = calcPe(data, curCell->neighCells[XP]->faces[XP]->u, curCell->faces[XP]->dy);
            f = data->rho * curCell->neighCells[XP]->faces[XP]->u;
            ae = D * curCell->faces[XP]->dy * A(Pe, data) + MAX(-f * curCell->faces[XP]->dy, 0.);
            // as
            D = data->alpha / (curCell->neighCells[XP]->x - curCell->x);
            Pe = calcPe(data, curCell->neighCells[YM]->faces[XP]->v, (curCell->neighCells[XP]->x - curCell->x));
            g = data->rho * curCell->neighCells[YM]->faces[XP]->v;
            as = D * (curCell->neighCells[XP]->x - curCell->x) * A(Pe, data) +
                MAX(g * (curCell->neighCells[XP]->x - curCell->x), 0.);
            // aw
            D = data->alpha / curCell->faces[XM]->dy;
            Pe = calcPe(data, curCell->faces[XM]->u, curCell->faces[XM]->dy);
            f = data->rho * curCell->faces[XM]->u;
            aw = D * curCell->faces[XM]->dy * A(Pe, data) + MAX(f * curCell->faces[XM]->dy, 0.);
            // ap
            ap = data->rho * curCell->faces[XM]->dy * (curCell->neighCells[XP]->x - curCell->x) / deltaT;
            // b
            b = ap * curCell->faces[XP]->u;
            // apTilde
            apTilde = an + ae + as + aw + ap;

            curCell->faces[XP]->u = an * curCell->neighCells[YP]->faces[XP]->u +
                ae * curCell->neighCells[XP]->faces[XP]->u + as * curCell->neighCells[YM]->faces[XP]->u +
                aw * curCell->faces[XM]->u + b;
            curCell->faces[XP]->u /= apTilde;

            // v
            // an
            D = data->alpha / curCell->faces[YP]->dx;
            Pe = calcPe(data, curCell->neighCells[YP]->faces[XP]->v, curCell->faces[YP]->dx);
            g = data->rho * curCell->neighCells[YP]->faces[XP]->v;
            an = D * curCell->faces[YP]->dx * A(Pe, data) + MAX(-g * curCell->faces[YP]->dx, 0.);
            // ae
            D = data->alpha / (curCell->neighCells[YP]->y - curCell->y);
            Pe = calcPe(data, curCell->neighCells[XP]->faces[XP]->u, (curCell->neighCells[YP]->y - curCell->y));
            f = data->rho * curCell->neighCells[XP]->faces[XP]->u;
            ae = D * (curCell->neighCells[YP]->y - curCell->y) * A(Pe, data) +
                MAX(-f * (curCell->neighCells[YP]->y - curCell->y), 0.);
            // as
            D = data->alpha / curCell->faces[YP]->dx;
            Pe = calcPe(data, curCell->neighCells[YM]->faces[XP]->v, curCell->faces[YP]->dx);
            g = data->rho * curCell->neighCells[YM]->faces[XP]->v;
            as = D * curCell->faces[YP]->dx * A(Pe, data) + MAX(g * curCell->faces[YP]->dx, 0.);
            // aw
            D = data->alpha / (curCell->neighCells[YP]->y - curCell->y);
            Pe = calcPe(data, curCell->faces[XM]->u, (curCell->neighCells[YP]->y - curCell->y));
            f = data->rho * curCell->faces[XM]->u;
            aw = D * (curCell->neighCells[YP]->y - curCell->y) * A(Pe, data) +
                MAX(f * (curCell->neighCells[YP]->y - curCell->y), 0.);
            // ap
            ap = data->rho * (curCell->neighCells[YP]->y - curCell->y) * curCell->faces[YP]->dx / deltaT;
            // b
            b = ap * curCell->faces[YP]->v;
            // apTilde
            apTilde = an + ae + as + aw + ap;

            curCell->faces[YP]->v = an * curCell->neighCells[YP]->faces[YP]->v +
                ae * curCell->neighCells[XP]->faces[YP]->v + as * curCell->faces[YM]->v +
                aw * curCell->neighCells[XM]->faces[YP]->v + b;
            curCell->faces[YP]->v /= apTilde;
        }
    }
}