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
bool solve(sData* data)
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

    while(curTime < data->maxTime && curIter < data->maxIter) {
        curIter++;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        if(data->solverType == CENTRAL)
            calcFluxCentral(data);
        else if(data->solverType == UPWIND)
            calcFluxUpwind(data);

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if(curCell->bType != 0)
                continue;
            ym = curCell->faces[YM];
            yp = curCell->faces[YP];
            xm = curCell->faces[XM];
            xp = curCell->faces[XP];

            // fluxes f* (horizontally)
            curCell->fluxBalance = ym->numFlux[M] * ym->dy - ym->numFlux[P] * ym->dx;
            curCell->fluxBalance += yp->numFlux[P] * yp->dx - yp->numFlux[M] * yp->dy;
            curCell->fluxBalance += xm->numFlux[M] * xm->dy - xm->numFlux[P] * xm->dx;
            curCell->fluxBalance += xp->numFlux[P] * xp->dx - xp->numFlux[M] * xp->dy;

            curCell->phi -= deltaT / (curCell->volume * data->rho) * curCell->fluxBalance;
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
            conv = data->rho * curFace->u + (curFace->neighCells[P]->phi + curFace->neighCells[M]->phi) / 2.;
            // diffusive part of f*
            dx = curFace->neighCells[P]->x - curFace->neighCells[M]->x;
            if(dx != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dx;
            // f* = f*(conv) + f*(diff)
            curFace->numFlux[M] = conv + diff;

            // convective part of g*
            conv = data->rho * curFace->v + (curFace->neighCells[P]->phi + curFace->neighCells[M]->phi) / 2.;
            // diffusive part of g*
            dy = curFace->neighCells[P]->y - curFace->neighCells[M]->y;
            if(dy != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dy;
            // g* = g*(conv) + g*(diff)
            curFace->numFlux[P] = conv + diff;
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
                conv = data->rho * curFace->u + curFace->neighCells[M]->phi;
            else
                conv = data->rho * curFace->u + curFace->neighCells[P]->phi;
            // diffusive part of f*
            dx = curFace->neighCells[P]->x - curFace->neighCells[M]->x;
            if(dx != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dx;
            // f* = f*(conv) + f*(diff)
            curFace->numFlux[M] = conv + diff;

            // convective part of g*
            if(curFace->u >= 0)
                conv = data->rho * curFace->v + curFace->neighCells[M]->phi;
            else
                conv = data->rho * curFace->v + curFace->neighCells[P]->phi;
            // diffusive part of g*
            dy = curFace->neighCells[P]->y - curFace->neighCells[M]->y;
            if(dy != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dy;
            // g* = g*(conv) + g*(diff)
            curFace->numFlux[P] = conv + diff;
        }
    }
}
