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
#include "setup.h"

//------------------------------------------------------
bool setup(sData* data)
{
    sCell* curCell = NULL;
    sFace* curFace = NULL;
    int i, j;

    /////////////////////////////////
    // construct mesh connectivity //
    /////////////////////////////////
    for(int cId = 0; cId < data->nCells; cId++) {
        // assign points to cells
        i = cId % (data->nPointsX - 1);
        j = cId / (data->nPointsX - 1);

        curCell = &data->cells[cId];
        curCell->points[XMYM] = &data->points[j * data->nPointsX + i];
        curCell->points[XPYM] = &data->points[j * data->nPointsX + i + 1];
        curCell->points[XMYP] = &data->points[(j + 1) * data->nPointsX + i];
        curCell->points[XPYP] = &data->points[(j + 1) * data->nPointsX + i + 1];

        if(cId == 0) {
            data->xMin = curCell->points[XMYM]->x;
            data->yMin = curCell->points[XMYM]->y;
        } else if(cId == data->nCells - 1) {
            data->xMax = curCell->points[XPYM]->x;
            data->yMax = curCell->points[XPYP]->y;
        }

        // assign faces to cells
        curCell->faces[YM] = &data->faces[cId];
        curCell->faces[YP] = &data->faces[cId + data->nPointsX - 1];
        curCell->faces[XM] = &data->faces[(data->nPointsX - 1) * data->nPointsY + cId + j];
        curCell->faces[XP] = &data->faces[(data->nPointsX - 1) * data->nPointsY + cId + j + 1];

        curCell->faces[YM]->id = cId;
        curCell->faces[YP]->id = cId + data->nCellsX;
        curCell->faces[XM]->id =
            cId + data->nCellsX * (data->nCellsY + 1) + (cId - (cId % data->nCellsX)) / data->nCellsX;
        curCell->faces[XP]->id =
            cId + data->nCellsX * (data->nCellsY + 1) + (cId - (cId % data->nCellsX)) / data->nCellsX + 1;

        // assign cells to faces
        curCell->faces[YM]->neighCells[P] = curCell;
        curCell->faces[YP]->neighCells[M] = curCell;
        curCell->faces[XM]->neighCells[P] = curCell;
        curCell->faces[XP]->neighCells[M] = curCell;

        // assign points to faces
        curCell->faces[YM]->points[M] = curCell->points[XMYM];
        curCell->faces[YM]->points[P] = curCell->points[XPYM];
        curCell->faces[YP]->points[M] = curCell->points[XMYP];
        curCell->faces[YP]->points[P] = curCell->points[XPYP];
        curCell->faces[XM]->points[M] = curCell->points[XMYM];
        curCell->faces[XM]->points[P] = curCell->points[XMYP];
        curCell->faces[XP]->points[M] = curCell->points[XPYM];
        curCell->faces[XP]->points[P] = curCell->points[XPYP];

        // assign neighboring cells to cells
        if(i != 0) {
            curCell->neighCells[XM] = &data->cells[cId - 1];
        }
        if(i != data->nPointsX - 2) {
            curCell->neighCells[XP] = &data->cells[cId + 1];
        }
        if(j != 0) {
            curCell->neighCells[YM] = &data->cells[cId - (data->nPointsX - 1)];
        }
        if(j != data->nPointsY - 2) {
            curCell->neighCells[YP] = &data->cells[cId + (data->nPointsX - 1)];
        }
    }

    /////////////////////////////////////
    // compute face centers and deltas //
    /////////////////////////////////////
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        curFace->dx = curFace->points[P]->x - curFace->points[M]->x;
        curFace->dy = curFace->points[P]->y - curFace->points[M]->y;
        curFace->x = (curFace->points[P]->x + curFace->points[M]->x) / 2.;
        curFace->y = (curFace->points[P]->y + curFace->points[M]->y) / 2.;
    }

    //////////////////////////////////////
    // compute cell centers and volumes //
    //////////////////////////////////////
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        curCell->x =
            curCell->points[XMYM]->x + curCell->points[XPYM]->x + curCell->points[XMYP]->x + curCell->points[XPYP]->x;
        curCell->x /= 4.;
        curCell->y =
            curCell->points[XMYM]->y + curCell->points[XPYM]->y + curCell->points[XMYP]->y + curCell->points[XPYP]->y;
        curCell->y /= 4.;

        curCell->volume = curCell->faces[XP]->dy * curCell->faces[YP]->dx;
    }

    ///////////////////////
    // set face velocity //
    ///////////////////////
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        curFace->u = data->u;
        curFace->v = data->v;
    }

    /////////////////////////////
    // set boundary conditions //
    /////////////////////////////
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];

        if((cId + data->nCellsX) % data->nCellsX == 0) { // linker Rand
            if(curCell->bTypeScalar == 1) {
                curCell->faces[XM]->bTypeScalar = 1;
                curCell->phi = curCell->bValueScalar;
            } else if(curCell->bTypeScalar == 2) {
                curCell->faces[XM]->bTypeScalar = 2;
                curCell->faces[XM]->bValueX = curCell->bValueScalarX;
                curCell->faces[XM]->bValueY = curCell->bValueScalarY;
            }
            if(curCell->bTypeVelocity == 1) {
                curCell->faces[XM]->bTypeVelocity = 1;
                curCell->faces[YM]->bTypeVelocity = 1;
                curCell->faces[YP]->bTypeVelocity = 1;
                curCell->faces[XM]->u = curCell->bValueU;
                curCell->faces[XM]->v = curCell->bValueV;
                curCell->faces[YM]->u = curCell->bValueU;
                curCell->faces[YM]->v = curCell->bValueV;
                curCell->faces[YP]->u = curCell->bValueU;
                curCell->faces[YP]->v = curCell->bValueV;
            } else if(curCell->bTypeVelocity == 2)
                curCell->faces[XM]->bTypeVelocity = 2;
            // Vpunkt = K * min(b,h)³ * max(b,h) / (12*eta*l) * deltaP mit K = 0.937 für b = 1, h = 1, l = 10
            // Vpunkt = 0.937 * 1³ * 1 / (12*1*10) * 10 = 0.937/12
            // u = Vpunkt/A = 0.937/12 / 1
            // curCell->faces[XM]->u = 0.937/12;
            // u = 1/(2*eta) * dpdx * (y² - y*h)
            double h = (data->yMax - data->yMin) * (1. - 1. / data->nCellsY);
            double ybar = curCell->y - (data->yMax - data->yMin) / (2. * data->nCellsY);
            double dpdx = (data->nCellsX - 1.) / (data->xMax - data->xMin);

            curCell->faces[XM]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
            curCell->faces[XM]->v = 0.;
            ybar = curCell->faces[YM]->y - (data->yMax - data->yMin) / (2 * data->nCellsY);
            curCell->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
            curCell->faces[YM]->v = 0.;
            ybar = curCell->faces[YP]->y - (data->yMax - data->yMin) / (2 * data->nCellsY);
            curCell->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
            curCell->faces[YP]->v = 0.;
        } else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) { // rechter Rand
            if(curCell->bTypeScalar == 1) {
                curCell->faces[XP]->bTypeScalar = 1;
                curCell->phi = curCell->bValueScalar;
            } else if(curCell->bTypeScalar == 2) {
                curCell->faces[XP]->bTypeScalar = 2;
                curCell->faces[XP]->bValueX = curCell->bValueScalarX;
                curCell->faces[XP]->bValueY = curCell->bValueScalarY;
            }
            if(curCell->bTypeVelocity == 1) {
                curCell->faces[XP]->bTypeVelocity = 1;
                curCell->faces[YM]->bTypeVelocity = 1;
                curCell->faces[YP]->bTypeVelocity = 1;
                curCell->faces[XP]->u = curCell->bValueU;
                curCell->faces[XP]->v = curCell->bValueV;
                curCell->faces[YM]->u = curCell->bValueU;
                curCell->faces[YM]->v = curCell->bValueV;
                curCell->faces[YP]->u = curCell->bValueU;
                curCell->faces[YP]->v = curCell->bValueV;
            } else if(curCell->bTypeVelocity == 2)
                curCell->faces[XP]->bTypeVelocity = 2;
            double h = (data->yMax - data->yMin) * (1. - 1. / data->nCellsY);
            double ybar = curCell->y - (data->yMax - data->yMin) / (2. * data->nCellsY);
            double dpdx = (data->nCellsX - 1.) / (data->xMax - data->xMin);

            curCell->faces[XP]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
            curCell->faces[XP]->v = 0.;
            ybar = curCell->faces[YM]->y - (data->yMax - data->yMin) / (2 * data->nCellsY);
            curCell->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
            curCell->faces[YM]->v = 0.;
            ybar = curCell->faces[YP]->y - (data->yMax - data->yMin) / (2 * data->nCellsY);
            curCell->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
            curCell->faces[YP]->v = 0.;
        }
        if(cId < data->nCellsX) { // unterer Rand
            if(curCell->bTypeScalar == 1) {
                curCell->faces[YM]->bTypeScalar = 1;
                curCell->phi = curCell->bValueScalar;
            } else if(curCell->bTypeScalar == 2) {
                curCell->faces[YM]->bTypeScalar = 2;
                curCell->faces[YM]->bValueX = curCell->bValueScalarX;
                curCell->faces[YM]->bValueY = curCell->bValueScalarY;
            }
            if(curCell->bTypeVelocity == 1) {
                curCell->faces[XM]->bTypeVelocity = 1;
                curCell->faces[XP]->bTypeVelocity = 1;
                curCell->faces[YM]->bTypeVelocity = 1;
                curCell->faces[XM]->u = curCell->bValueU;
                curCell->faces[XM]->v = curCell->bValueV;
                curCell->faces[XP]->u = curCell->bValueU;
                curCell->faces[XP]->v = curCell->bValueV;
                curCell->faces[YM]->u = curCell->bValueU;
                curCell->faces[YM]->v = curCell->bValueV;
            } else if(curCell->bTypeVelocity == 2)
                curCell->faces[YM]->bTypeVelocity = 2;
        } else if(data->nCells - (cId + 1) < data->nCellsX) { // oberer Rand
            if(curCell->bTypeScalar == 1) {
                curCell->faces[YP]->bTypeScalar = 1;
                curCell->phi = curCell->bValueScalar;
            } else if(curCell->bTypeScalar == 2) {
                curCell->faces[YP]->bTypeScalar = 2;
                curCell->faces[YP]->bValueX = curCell->bValueScalarX;
                curCell->faces[YP]->bValueY = curCell->bValueScalarY;
            }
            if(curCell->bTypeVelocity == 1) {
                curCell->faces[XM]->bTypeVelocity = 1;
                curCell->faces[XP]->bTypeVelocity = 1;
                curCell->faces[YP]->bTypeVelocity = 1;
                curCell->faces[XM]->u = curCell->bValueU;
                curCell->faces[XM]->v = curCell->bValueV;
                curCell->faces[XP]->u = curCell->bValueU;
                curCell->faces[XP]->v = curCell->bValueV;
                curCell->faces[YP]->u = curCell->bValueU;
                curCell->faces[YP]->v = curCell->bValueV;
            }
        } else if(curCell->bTypeVelocity == 2)
            curCell->faces[YP]->bTypeVelocity = 2;
        double h = (data->yMax - data->yMin) * (1. - 1. / data->nCellsY);
        double ybar = curCell->y - (data->yMax - data->yMin) / (2. * data->nCellsY);
        double dpdx = (data->nCellsX - 1.) / (data->xMax - data->xMin);

        curCell->faces[XM]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
        curCell->faces[XM]->v = 0.;
        ybar = curCell->faces[YM]->y - (data->yMax - data->yMin) / (2 * data->nCellsY);
        curCell->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
        curCell->faces[YM]->v = 0.;
        ybar = curCell->faces[YP]->y - (data->yMax - data->yMin) / (2 * data->nCellsY);
        curCell->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (ybar * ybar - ybar * h);
        curCell->faces[YP]->v = 0.;
    }

    return true;
}

void setRigidBodyBoundaries(sData* data)
{

    double omega = 5.;
    double rSquared;

    sCell* curCell;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        if(curCell->bTypeVelocity == 1) {
            rSquared = curCell->x * curCell->x + curCell->y * curCell->y;
            curCell->p = 100. + 31.25 * rSquared;
            // data->cells[cId].p = 100. + 0. * ((cId+data->nCellsX) % (data->nCellsX));
        }
    }

    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        if(curFace->bTypeVelocity == 1) {
            curFace->u = -curFace->y * omega;
            curFace->v = +curFace->x * omega;
        }
    }
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        std::cout << fId << ":   " << curFace->u << "   " << curFace->v << std::endl;
    }
    std::cout << "\nPressure:" << std::endl;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        std::cout << cId << ":   " << curCell->p << std::endl;
    }
}