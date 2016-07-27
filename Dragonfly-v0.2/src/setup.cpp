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
            } else if(curCell->bTypeVelocity == 2) {
                curCell->faces[XM]->bTypeVelocity = 2;
                curCell->faces[YM]->bTypeVelocity = 2;
                curCell->faces[YP]->bTypeVelocity = 2;
            }
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
            } else if(curCell->bTypeVelocity == 2) {
                curCell->faces[XP]->bTypeVelocity = 2;
                curCell->faces[YM]->bTypeVelocity = 2;
                curCell->faces[YP]->bTypeVelocity = 2;
            }
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
            } else if(curCell->bTypeVelocity == 2) {
                curCell->faces[XM]->bTypeVelocity = 2;
                curCell->faces[YM]->bTypeVelocity = 2;
                curCell->faces[XP]->bTypeVelocity = 2;
            }
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
            } else if(curCell->bTypeVelocity == 2) {
                curCell->faces[XM]->bTypeVelocity = 2;
                curCell->faces[YP]->bTypeVelocity = 2;
                curCell->faces[XP]->bTypeVelocity = 2;
            }
        }
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
        
        if((cId + data->nCellsX) % data->nCellsX == 0) // linker Rand
            curCell->bTypePressure = DIRICHLET;
        else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) // rechter Rand
            curCell->bTypePressure = DIRICHLET;
        else if(data->nCells - (cId + 1) < data->nCellsX) // oberer Rand
            curCell->bTypePressure = DIRICHLET;
        else if(cId < data->nCellsX) // unterer Rand
            curCell->bTypePressure = DIRICHLET; 
            
        if(curCell->bTypePressure == DIRICHLET) {
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
}

void setCouetteBoundaries(sData* data)
{
    double U = 5.;
    double Uin = 2.5;
    double pressure = 10.;

    sCell* curCell;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        
        if((cId + data->nCellsX) % data->nCellsX == 0) // linker Rand
            curCell->bTypePressure = NEUMANN;
        else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) // rechter Rand
            curCell->bTypePressure = DIRICHLET;
        else if(data->nCells - (cId + 1) < data->nCellsX) // oberer Rand
            curCell->bTypePressure = NEUMANN;
        else if(cId < data->nCellsX) // unterer Rand
            curCell->bTypePressure = NEUMANN; 
        else
           curCell->bTypePressure = INNERCELL;
        
        curCell->p = pressure;
        if (curCell->bTypePressure == INNERCELL) {
           curCell->p = 1000.;
        }
    }

    int nX = data->nCellsX;
    int nY = data->nCellsY;

    double b = data->yMax-data->yMin;
    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
         /*curFace->u = U/b*curFace->y;
         curFace->v = 0.;
         continue;*/

        if(fId == 0) { // UL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX - 1) { // UR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * nY) { // OL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[XP]->u = U;
        } else if(fId == nX * (nY + 1) - 1) { // OR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[XM]->u = U;
        } else if(fId == nX * (nY + 1)) { // LU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            /*curFace->u = U/b*curFace->y;
            curFace->neighCells[P]->faces[YP]->u = U/b*curFace->neighCells[P]->faces[YP]->y;*/
            curFace->u = Uin;
            curFace->neighCells[P]->faces[YP]->u = Uin;
        } else if(fId == nX * (nY + 2)) { // RU
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = NEUMANN;
            /*curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[YP]->u = Uin;*/
        } else if(fId == (2 * nX + 1) * nY - 1) { // LO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->u = Uin;
            curFace->neighCells[P]->faces[YM]->u = Uin;
            /*curFace->u = U/b*curFace->y;
            curFace->neighCells[P]->faces[YM]->u = U/b*curFace->neighCells[P]->faces[YM]->y;*/
        } else if(fId == (2 * nX + 1) * nY - 1 + nX) { // RO
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = NEUMANN;
            /*curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[YM]->u = Uin;*/
        } else if(fId <= nX - 1) { // U
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if((fId - (nX * (nY + 1))) % (nX + 1) == 0 && (fId - (nX * (nY + 1)))>=0 ) { // L
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            curFace->u = Uin;
            curFace->neighCells[P]->faces[YM]->u = Uin;
            curFace->neighCells[P]->faces[YP]->u = Uin;
            /*curFace->u = U/b*curFace->y;
            curFace->neighCells[P]->faces[YM]->u = U/b*curFace->neighCells[P]->faces[YM]->y;
            curFace->neighCells[P]->faces[YP]->u = U/b*curFace->neighCells[P]->faces[YP]->y;*/
        } else if((fId - (nX * (nY + 2))) % (nX + 1) == 0 && (fId - (nX * (nY + 2)))>=0 ) { // R
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = NEUMANN;
            /*curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
            curFace->u = Uin;
            curFace->neighCells[M]->faces[YM]->u = Uin;
            curFace->neighCells[M]->faces[YP]->u = Uin;*/
        } else if(fId >= nX * nY && fId <= nX * (nY + 1) - 1) { // O
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[XM]->u = U;
            curFace->neighCells[M]->faces[XP]->u = U;
        }
    }
}

void setPoiseuilleBoundaries(sData* data)
{
    int nX = data->nCellsX;
    int nY = data->nCellsY;

    sCell* curCell;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        curCell->p = (99 - cId % nX) * 100. / 99.;
    }

    sFace* curFace;
    double dpdx = (data->nCellsX - 1.) / (data->xMax - data->xMin);
    double b = data->yMax - data->yMin;
    double y;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(fId == 0) { // UL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX - 1) { // UR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * nY) { // OL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 1) - 1) { // OR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 1)) { // LU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YP]->y;
            curFace->neighCells[P]->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId == nX * (nY + 2)) { // RU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[M]->faces[YP]->y;
            curFace->neighCells[M]->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId == (2 * nX + 1) * nY - 1) { // LO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YM]->y;
            curFace->neighCells[P]->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId == (2 * nX + 1) * nY - 1 + nX) { // RO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[M]->faces[YM]->y;
            curFace->neighCells[M]->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId <= nX - 1) { // U
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
            y = curFace->neighCells[P]->faces[XM]->y;
            curFace->neighCells[P]->faces[XM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[XP]->y;
            curFace->neighCells[P]->faces[XP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if((fId - (nX * (nY + 1))) % (nX + 1) == 0 && (fId - (nX * (nY + 1))) >= 0) { // L
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YM]->y;
            curFace->neighCells[P]->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YP]->y;
            curFace->neighCells[P]->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if((fId - (nX * (nY + 2))) % (nX + 1) == 0 && (fId - (nX * (nY + 2))) >= 0) { // R
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[M]->faces[YM]->y;
            curFace->neighCells[M]->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[M]->faces[YP]->y;
            curFace->neighCells[M]->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId >= nX * nY && fId <= nX * (nY + 1) - 1) { // O
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
            y = curFace->neighCells[M]->faces[XM]->y;
            curFace->neighCells[M]->faces[XM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[M]->faces[XP]->y;
            curFace->neighCells[M]->faces[XP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else {
            curFace->u = 0.;
            curFace->v = 0.;
        }
    }

    for(int cId = 0; cId < data->nCells; ++cId) {
        std::cout << cId << ": " << data->cells[cId].p << std::endl;
    }
    std::cout << std::endl;

    for(int fId = 0; fId < data->nFaces; ++fId) {
        std::cout << fId << ": " << data->faces[fId].u << std::endl;
    }
    std::cout << std::endl;
}

void setDrivenCavityBoundaries(sData* data)
{
    double U = 9./3.;
    double pressure = 100.;
    
    int nX = data->nCellsX;
    int nY = data->nCellsY;

    sCell* curCell;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        
        curCell->p = pressure;
        
        if((cId + data->nCellsX) % data->nCellsX == 0) // linker Rand
            curCell->bTypePressure = NEUMANN;
        else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) // rechter Rand
            curCell->bTypePressure = NEUMANN;
        if(data->nCells - (cId + 1) < data->nCellsX) { // oberer Rand
            curCell->bTypePressure = DIRICHLET;
            curCell->p = pressure;
        }
        else if(cId < data->nCellsX) // unterer Rand
            curCell->bTypePressure = NEUMANN;    
    }

    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(fId == 0) { // UL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX - 1) { // UR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * nY) { // OL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[XP]->u = U;
        } else if(fId == nX * (nY + 1) - 1) { // OR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[XM]->u = U;
        } else if(fId == nX * (nY + 1)) { // LU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 2)) { // RU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
        } else if(fId == (2 * nX + 1) * nY - 1) { // LO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
        } else if(fId == (2 * nX + 1) * nY - 1 + nX) { // RO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
        } else if(fId <= nX - 1) { // U
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if((fId - (nX * (nY + 1))) % (nX + 1) == 0 && (fId - (nX * (nY + 1))) >= 0) { // L
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
        } else if((fId - (nX * (nY + 2))) % (nX + 1) == 0 && (fId - (nX * (nY + 2))) >= 0) { // R
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
        } else if(fId >= nX * nY && fId <= nX * (nY + 1) - 1) { // O
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[M]->faces[XM]->u = U;
            curFace->neighCells[M]->faces[XP]->u = U;
        } else {
            curFace->u = 0.;
            curFace->v = 0.;
        }
    }
}

void setRiverFlowBoundariers1(sData * data)
{
    double U = 5.;
    double V = 3.;
    double source = 100.;

    int tmp = (int)(data->nCellsX / 3 + 1);
    int iter = 0;
    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        curFace->u = U;
        curFace->v = 0.;
        if((fId + tmp) % (data->nCellsX + 1) == 0 && curFace->dx == 0 && curFace->y < 0.3) {
            curFace->v = V;
            curFace->neighCells[M]->faces[YM]->v = V;
            if(iter == 1)
                curFace->neighCells[M]->s = source;
            ++iter;
        }
    }
}

void setRiverFlowBoundariers2(sData * data)
{
    double U = 5.;
    double source1 = 10.;
    double source2 = -100.;

    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        curFace->u = U;
        curFace->v = 0.;
        if(fId == 7540)
            curFace->neighCells[P]->s = source1;
        else if(fId == 7550)
            curFace->neighCells[M]->s = source2;
    }
}

void setRiverFlowBoundariers3(sData * data)
{
    double U = 5.;
    double Vin = 3.;
    double Vout = -5.;
    double source1 = 100.;
    double source2 = -200.;

    int iter1 = 0;
    int iter2 = 0;
    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];
        curFace->u = U;
        curFace->v = 0.;
        if((fId - 7540) % (data->nCellsX + 1) == 0 && curFace->y < 0.3 && curFace->dx == 0) {
            curFace->v = Vin;
            curFace->neighCells[M]->faces[YM]->v = Vin;
            if(iter1 == 1) {
                curFace->neighCells[M]->s = source1;
            }
            ++iter1;
        }
        if((fId - 7550) % (data->nCellsX + 1) == 0 && curFace->y < 0.5 && curFace->dx == 0) {
            curFace->v = Vout;
            curFace->neighCells[M]->faces[YM]->v = Vout;
            if(iter2 == 1) {
                curFace->neighCells[M]->s = source2;
            }
            ++iter2;
        }
    }
}

void setGartenschlauchBoundaries(sData * data)
{
    double U = 5.;
    double pressure = 100.;
    
    int nX = data->nCellsX;
    int nY = data->nCellsY;

    sCell* curCell;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
                
        if((cId + data->nCellsX) % data->nCellsX == 0) // linker Rand
            curCell->bTypePressure = NEUMANN;
        else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) { // rechter Rand
            curCell->bTypePressure = DIRICHLET;
            curCell->p = pressure;
        }
        else if(data->nCells - (cId + 1) < data->nCellsX)  // oberer Rand
            curCell->bTypePressure = NEUMANN;
        else if(cId < data->nCellsX) // unterer Rand
            curCell->bTypePressure = NEUMANN;    
    }

    sFace* curFace;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(fId == 0) { // UL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX - 1) { // UR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * nY) { // OL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 1) - 1) { // OR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 1)) { // LU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[P]->faces[YP]->u  = U;
        } else if(fId == nX * (nY + 2)) { // RU
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = NEUMANN;
        } else if(fId == (2 * nX + 1) * nY - 1) { // LO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[P]->faces[YM]->u  = U;
        } else if(fId == (2 * nX + 1) * nY - 1 + nX) { // RO
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = NEUMANN;
        } else if(fId <= nX - 1) { // U
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if((fId - (nX * (nY + 1))) % (nX + 1) == 0 && (fId - (nX * (nY + 1))) >= 0) { // L
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            curFace->u = U;
            curFace->neighCells[P]->faces[YM]->u  = U;
            curFace->neighCells[P]->faces[YP]->u  = U;
        } else if((fId - (nX * (nY + 2))) % (nX + 1) == 0 && (fId - (nX * (nY + 2))) >= 0) { // R
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = NEUMANN;
        } else if(fId >= nX * nY && fId <= nX * (nY + 1) - 1) { // O
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else {
            curFace->u = 0.;
            curFace->v = 0.;
        }
    }
}

void setStolperdrahtBoundaries(sData* data) {
    int nX = data->nCellsX;
    int nY = data->nCellsY;
    double pressure = 10.;

    sCell* curCell;
    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
                
        if((cId + data->nCellsX) % data->nCellsX == 0) // linker Rand
            curCell->bTypePressure = NEUMANN;
        else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) { // rechter Rand
            curCell->bTypePressure = DIRICHLET;
            curCell->p = pressure;
        }
        else if(data->nCells - (cId + 1) < data->nCellsX)  // oberer Rand
            curCell->bTypePressure = NEUMANN;
        else if(cId < data->nCellsX) // unterer Rand
            curCell->bTypePressure = NEUMANN;   

         if (cId==20 || cId==21 || cId==70 || cId==71) {
            curCell->bTypeVelocity = DIRICHLET;
            curCell->faces[XM]->bTypeVelocity = DIRICHLET;
            curCell->faces[XP]->bTypeVelocity = DIRICHLET;
            curCell->faces[YP]->bTypeVelocity = DIRICHLET;
            curCell->faces[XM]->u = 0.;
            curCell->faces[XM]->v = 0.;
            curCell->faces[YP]->u = 0.;
            curCell->faces[YP]->v = 0.;
            curCell->faces[XP]->u = 0.;
            curCell->faces[XP]->v = 0.;
            curCell->bTypePressure = DIRICHLET;
         }
    }

    sFace* curFace;
    double dpdx = 10.;
    double b = data->yMax - data->yMin;
    double y;
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(fId == 0) { // UL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX - 1) { // UR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * nY) { // OL
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 1) - 1) { // OR
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
        } else if(fId == nX * (nY + 1)) { // LU
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YP]->y;
            curFace->neighCells[P]->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId == nX * (nY + 2)) { // RU
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = NEUMANN;
        } else if(fId == (2 * nX + 1) * nY - 1) { // LO
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YM]->y;
            curFace->neighCells[P]->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if(fId == (2 * nX + 1) * nY - 1 + nX) { // RO
            curFace->bTypeVelocity = NEUMANN;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = NEUMANN;
        } else if(fId <= nX - 1) { // U
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else if((fId - (nX * (nY + 1))) % (nX + 1) == 0 && (fId - (nX * (nY + 1))) >= 0) { // L
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[P]->faces[YP]->bTypeVelocity = DIRICHLET;
            y = curFace->y;
            curFace->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YM]->y;
            curFace->neighCells[P]->faces[YM]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
            y = curFace->neighCells[P]->faces[YP]->y;
            curFace->neighCells[P]->faces[YP]->u = -1. / (2. * data->eta) * dpdx * (y * y - y * b);
        } else if((fId - (nX * (nY + 2))) % (nX + 1) == 0 && (fId - (nX * (nY + 2))) >= 0) { // R
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[YP]->bTypeVelocity = DIRICHLET;
        } else if(fId >= nX * nY && fId <= nX * (nY + 1) - 1) { // O
            curFace->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XM]->bTypeVelocity = DIRICHLET;
            curFace->neighCells[M]->faces[XP]->bTypeVelocity = DIRICHLET;
        } else {
            curFace->u = 0.;
            curFace->v = 0.;
        }
    }   
}