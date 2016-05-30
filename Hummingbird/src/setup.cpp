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

//-----------------------------------------------------
bool setup(sData* data)
{
    std::cout << "\nSetup:\n-------\n";

    // setup xi/et grid
    for(int i = 0; i < data->nX; ++i)
        for(int j = 0; j < data->nY; ++j) {
            data->xi[i][j] = i;
            data->et[i][j] = j;
        }
    printValue(data, "xi", data->xi);
    printValue(data, "et", data->et);

    // setup x/y grid
    for(int i = 0; i < data->nX; ++i)
        for(int j = 0; j < data->nY; ++j) {
            data->x[i][j] = X(data, i, j);
            data->y[i][j] = Y(data, i, j);
        }
    printValue(data, "x", data->x);
    printValue(data, "y", data->y);

    // set inital values of scalars s1 and s2
    for (int i = 0; i < data->nX; ++i) {
       for (int j = 0; j < data->nY; ++j) {
          data->s1[i][j] = data->uInfty*data->x[i][j]+data->vInfty*data->y[i][j];
          data->s2[i][j] = -data->vInfty*data->x[i][j]+data->uInfty*data->y[i][j];
       }
    }

    // precompute derivatives for coordinate transformation
    double x, y, dx, dy;
    for(int i = 1; i < data->nX - 1; i++) {
        for(int j = 1; j < data->nY - 1; j++) {
            // calculate derivatives at inner grid points
            x = data->x[i][j];
            y = data->y[i][j];
            dx = 1. / DELD;
            dy = 1. / DELD;
            data->dxidx[i][j] = (Xi(data, x + dx, y) - Xi(data, x - dx, y)) / (2 * dx);
            data->detdx[i][j] = (Et(data, x + dx, y) - Et(data, x - dx, y)) / (2 * dx);
            data->detdy[i][j] = (Et(data, x, y + dy) - Et(data, x, y - dy)) / (2 * dy);
            dx = 1. / DELDD;
            dy = 1. / DELDD;
            data->ddxidx[i][j] = (Xi(data, x + dx, y) + Xi(data, x - dx, y) - 2 * Xi(data, x, y)) / (dx * dx);
            data->ddetdx[i][j] = (Et(data, x + dx, y) + Et(data, x - dx, y) - 2 * Et(data, x, y)) / (dx * dx);
            data->ddetdy[i][j] = (Et(data, x, y + dy) + Et(data, x, y - dy) - 2 * Et(data, x, y)) / (dy * dy);
        }
        // at lower boundary
        x = data->x[i][0];
        y = data->y[i][0];
        dx = 1. / DELD;
        data->dxidx[i][0] = (Xi(data, x + dx, y) - Xi(data, x - dx, y)) / (2 * dx);
        data->detdx[i][0] = (Et(data, x + dx, y) - Et(data, x - dx, y)) / (2 * dx);
        // at  upper boundary
        x = data->x[i][data->nY - 1];
        y = data->y[i][data->nY - 1];
        dx = 1. / DELD;
        data->dxidx[i][data->nY - 1] = (Xi(data, x + dx, y) - Xi(data, x - dx, y)) / (2 * dx);
        data->detdx[i][data->nY - 1] = (Et(data, x + dx, y) - Et(data, x - dx, y)) / (2 * dx);
    }
    for(int j = 1; j < data->nY - 1; ++j) {
        // at left boundary
        x = data->x[0][j];
        y = data->y[0][j];
        dy = 1. / DELD;
        data->dxidy[0][j] = (Xi(data, x, y + dy) - Xi(data, x, y - dy)) / (2 * dy);
        data->detdy[0][j] = (Et(data, x, y + dy) - Et(data, x, y - dy)) / (2 * dy);
        // at  right boundary
        x = data->x[data->nX - 1][j];
        y = data->y[data->nX - 1][j];
        dy = 1. / DELD;
        data->detdy[data->nX - 1][j] = (Et(data, x, y + dy) - Et(data, x, y - dy)) / (2 * dy);
    }

    for(int j = 0; j < data->nY; ++j) {
        // at left boundary
        x = data->x[0][j];
        y = data->y[0][j];
        dx = 1. / DELD;
        data->dxidx[0][j] = (Xi(data, x + dx, y) - Xi(data, x, y)) / dx;
        data->detdx[0][j] = (Et(data, x + dx, y) - Et(data, x, y)) / dx;
        // at  right boundary
        x = data->x[data->nX - 1][j];
        y = data->y[data->nX - 1][j];
        dx = 1. / DELD;
        data->dxidx[data->nX - 1][j] = (Xi(data, x, y) - Xi(data, x - dx, y)) / dx;
        data->detdx[data->nX - 1][j] = (Et(data, x, y) - Et(data, x - dx, y)) / dx;
    }

    for(int i = 0; i < data->nX; ++i) {
        // at lower boundary
        x = data->x[i][0];
        y = data->y[i][0];
        dy = 1. / DELD;
        data->dxidy[i][0] = (Xi(data, x, y + dy) - Xi(data, x, y)) / dy;
        data->detdy[i][0] = (Et(data, x, y + dy) - Et(data, x, y)) / dy;
        // at upper boundary
        x = data->x[i][data->nY - 1];
        y = data->y[i][data->nY - 1];
        dy = 1. / DELD;
        data->dxidy[i][data->nY - 1] = (Xi(data, x, y) - Xi(data, x, y - dy)) / dy;
        data->detdy[i][data->nY - 1] = (Et(data, x, y) - Et(data, x, y - dy)) / dy;
    }

    printValue(data, "dxidx", data->dxidx);
    printValue(data, "detdx", data->detdx);
    printValue(data, "dxidy", data->dxidy);
    printValue(data, "detdy", data->detdy);
    printValue(data, "ddxidx", data->ddxidx);
    printValue(data, "ddetdx", data->ddetdx);
    printValue(data, "ddxidy", data->ddxidy);
    printValue(data, "ddetdy", data->ddetdy);

    for(int i = 1; i < data->nX - 1; i++) {
        for(int j = 1; j < data->nY - 1; j++) {

            data->a1[i][j] = data->dxidx[i][j] * data->dxidx[i][j] + data->dxidy[i][j] * data->dxidy[i][j];
            data->a2[i][j] = data->detdx[i][j] * data->detdx[i][j] + data->detdy[i][j] * data->detdy[i][j];
            data->a3[i][j] = 2 * (data->dxidx[i][j] * data->detdx[i][j] + data->dxidy[i][j] * data->detdy[i][j]);
            data->a4[i][j] = data->ddxidx[i][j] + data->ddxidy[i][j];
            data->a5[i][j] = data->ddetdx[i][j] + data->ddetdy[i][j];
        }
    }

    // set boundary conditions for scalar fields
    if(data->potentialFunc == CHANNELFLOW) {
        double dhdx[2];

        data->uOut = (ContourMax(data->xMin) - ContourMin(data->xMin)) /
            (ContourMax(data->xMax) - ContourMin(data->xMax)) * data->uIn;

        for(int j = 0; j < data->nY; ++j) {
            // left boundary
            data->s1[0][j] = data->uIn * data->x[0][j];
            data->s2[0][j] = data->uIn * data->y[0][j];

            // right boundary
            data->s1[data->nX - 1][j] = data->uOut * data->x[data->nX - 1][j];
            // data->s2[data->nX - 1][j] = data->uOut * data->y[data->nX - 1][j];
        }
        for(int j = 0; j < data->nY; ++j) {
            data->s2[data->nX - 1][j] = data->s2[0][0] +
                (data->s2[0][data->nY - 1] - data->s2[0][0]) /
                    (data->y[data->nX - 1][data->nY - 1] - data->y[data->nX - 1][0]) *
                    (data->y[data->nX - 1][j] - data->y[data->nX - 1][0]);
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
                dhdx[1] *
                ((data->s1[i + 1][data->nY - 1] - data->s1[i - 1][data->nY - 1]) * data->dxidx[i][data->nY - 1] +
                    (data->s1[i][data->nY - 3] - 4 * data->s1[i][data->nY - 2]) * data->detdx[i][data->nY - 1]);
            data->s1[i][data->nY - 1] -=
                (data->s1[i + 1][data->nY - 1] - data->s1[i - 1][data->nY - 1]) * data->dxidy[i][data->nY - 1] +
                (data->s1[i][data->nY - 3] - 4 * data->s1[i][data->nY - 2]) * data->detdy[i][data->nY - 1];
            data->s1[i][data->nY - 1] /= 3 * (data->detdy[i][data->nY - 1] - data->detdx[i][data->nY - 1] * dhdx[1]);
            data->s2[i][data->nY - 1] = data->s2[0][data->nY - 1];
        }
    } else if(data->potentialFunc == PARALLELFLOW) {
        for(int i = 0; i < data->nX; ++i) {
            data->s1[i][0] = data->uInfty * data->x[i][0] + data->vInfty * data->y[i][0];
            data->s2[i][0] = -data->vInfty * data->x[i][0] + data->uInfty * data->y[i][0];
            data->s1[i][data->nY - 1] =
                data->uInfty * data->x[i][data->nY - 1] + data->vInfty * data->y[i][data->nY - 1];
            data->s2[i][data->nY - 1] =
                -data->vInfty * data->x[i][data->nY - 1] + data->uInfty * data->y[i][data->nY - 1];
        }
        for(int j = 0; j < data->nY; ++j) {
            data->s1[0][j] = data->uInfty * data->x[0][j] + data->vInfty * data->y[0][j];
            data->s2[0][j] = -data->vInfty * data->x[0][j] + data->uInfty * data->y[0][j];
            data->s1[data->nX - 1][j] =
                data->uInfty * data->x[data->nX - 1][j] + data->vInfty * data->y[data->nX - 1][j];
            data->s2[data->nX - 1][j] =
                -data->vInfty * data->x[data->nX - 1][j] + data->uInfty * data->y[data->nX - 1][j];
        }
    } else if(data->potentialFunc == STAGNATION_POINT) {
        for(int i = 0; i < data->nX; ++i) {
            data->s1[i][0] = data->magnitude * (data->x[i][0] * data->x[i][0] - data->y[i][0] * data->y[i][0]);
            data->s2[i][0] = 2. * data->magnitude * data->x[i][0] * data->y[i][0];
            data->s1[i][data->nY - 1] = data->magnitude * (data->x[i][data->nY - 1] * data->x[i][data->nY - 1] -
                                                              data->y[i][data->nY - 1] * data->y[i][data->nY - 1]);
            data->s2[i][data->nY - 1] = 2. * data->magnitude * data->x[i][data->nY - 1] * data->y[i][data->nY - 1];
        }
        for(int j = 0; j < data->nY; ++j) {
            data->s1[0][j] = data->magnitude * (data->x[0][j] * data->x[0][j] - data->y[0][j] * data->y[0][j]);
            data->s2[0][j] = 2. * data->magnitude * data->x[0][j] * data->y[0][j];
            data->s1[data->nX - 1][j] = data->magnitude * (data->x[data->nX - 1][j] * data->x[data->nX - 1][j] -
                                                              data->y[data->nX - 1][j] * data->y[data->nX - 1][j]);
            data->s2[data->nX - 1][j] = 2. * data->magnitude * data->x[data->nX - 1][j] * data->y[data->nX - 1][j];
        }
    } else if(data->potentialFunc == SOURCE) {
        for(int i = 0; i < data->nX; ++i) {
            data->s1[i][0] = data->magnitude / (2 * PI) *
                logf(sqrt((data->x[i][0] * data->x[i][0]) + (data->y[i][0] * data->y[i][0])));
            data->s1[i][data->nY - 1] =
                data->magnitude / (2 * PI) * logf(sqrt((data->x[i][data->nY - 1] * data->x[i][data->nY - 1]) +
                                                 (data->y[i][data->nY - 1] * data->y[i][data->nY - 1])));
            data->s2[i][0] = data->magnitude / (2 * PI) * atanf(data->y[i][0] / data->x[i][0]);
            data->s2[i][data->nY - 1] =
                data->magnitude / (2 * PI) * atanf(data->y[i][data->nY - 1] / data->x[i][data->nY - 1]);
        }
        for(int j = 0; j < data->nY; ++j) {
            data->s1[0][j] = data->magnitude / (2 * PI) *
                logf(sqrt((data->x[0][j] * data->x[0][j]) + (data->y[0][j] * data->y[0][j])));
            data->s1[data->nX - 1][j] =
                data->magnitude / (2 * PI) * logf(sqrt((data->x[data->nX - 1][j] * data->x[data->nX - 1][j]) +
                                                 (data->y[data->nX - 1][j] * data->y[data->nX - 1][j])));
            data->s2[0][j] = data->magnitude / (2 * PI) * atanf(data->y[0][j] / data->x[0][j]);
            data->s2[data->nX - 1][j] =
                data->magnitude / (2 * PI) * atanf(data->y[data->nX - 1][j] / data->x[data->nX - 1][j]);
        }
    } else if(data->potentialFunc == POTENTIAL_VORTEX) {
        for(int i = 0; i < data->nX; ++i) {
            data->s1[i][0] = data->magnitude / (2 * PI) * atanf(data->y[i][0] / data->x[i][0]);
            data->s1[i][data->nY - 1] =
                data->magnitude / (2 * PI) * atanf(data->y[i][data->nY - 1] / data->x[i][data->nY - 1]);
            data->s2[i][0] = -data->magnitude / (2 * PI) *
                logf(sqrt((data->x[i][0] * data->x[i][0]) + (data->y[i][0] * data->y[i][0])));
            data->s2[i][data->nY - 1] =
                -data->magnitude / (2 * PI) * logf(sqrt((data->x[i][data->nY - 1] * data->x[i][data->nY - 1]) +
                                                 (data->y[i][data->nY - 1] * data->y[i][data->nY - 1])));
        }
        for(int j = 0; j < data->nY; ++j) {
            data->s1[0][j] = data->magnitude / (2 * PI) * atanf(data->y[0][j] / data->x[0][j]);
            data->s1[data->nX - 1][j] =
                data->magnitude / (2 * PI) * atanf(data->y[data->nX - 1][j] / data->x[data->nX - 1][j]);
            data->s2[0][j] = -data->magnitude / (2 * PI) *
                logf(sqrt((data->x[0][j] * data->x[0][j]) + (data->y[0][j] * data->y[0][j])));
            data->s2[data->nX - 1][j] =
                -data->magnitude / (2 * PI) * logf(sqrt((data->x[data->nX - 1][j] * data->x[data->nX - 1][j]) +
                                                 (data->y[data->nX - 1][j] * data->y[data->nX - 1][j])));
        }
    } else if(data->potentialFunc == DIPOL) {
        for(int i = 0; i < data->nX; ++i) {
            data->s1[i][0] = data->magnitude / (2 * PI) * data->x[i][0] /
                sqrtf(data->x[i][0] * data->x[i][0] + data->y[i][0] * data->y[i][0]);
            data->s1[i][data->nY - 1] = data->magnitude / (2 * PI) * data->x[i][data->nY - 1] /
                sqrtf(data->x[i][data->nY - 1] * data->x[i][data->nY - 1] +
                                            data->y[i][data->nY - 1] * data->y[i][data->nY - 1]);
            data->s2[i][0] = -data->magnitude / (2 * PI) * data->y[i][0] /
                sqrtf(data->x[i][0] * data->x[i][0] + data->y[i][0] * data->y[i][0]);
            data->s2[i][data->nY - 1] = -data->magnitude / (2 * PI) * data->y[i][data->nY - 1] /
                sqrtf(data->x[i][data->nY - 1] * data->x[i][data->nY - 1] +
                                            data->y[i][data->nY - 1] * data->y[i][data->nY - 1]);
        }
        for(int j = 0; j < data->nY; ++j) {
            data->s1[0][j] = data->magnitude / (2 * PI) * data->x[0][j] /
                sqrtf(data->x[0][j] * data->x[0][j] + data->y[0][j] * data->y[0][j]);
            data->s1[data->nX - 1][j] = data->magnitude / (2 * PI) * data->x[data->nX - 1][j] /
                sqrtf(data->x[data->nX - 1][j] * data->x[data->nX - 1][j] + data->y[data->nX - 1][j] * data->y[data->nX - 1][j]);
            data->s2[0][j] = -data->magnitude / (2 * PI) * data->y[0][j] /
                sqrtf(data->x[0][j] * data->x[0][j] + data->y[0][j] * data->y[0][j]);
            data->s2[data->nX - 1][j] = -data->magnitude / (2 * PI) * data->y[data->nX - 1][j] /
                sqrtf(data->x[data->nX - 1][j] * data->x[data->nX - 1][j] + data->y[data->nX - 1][j] * data->y[data->nX - 1][j]);
        }
    }

    printValue(data, "s boundaries", data->s1);
    printValue(data, "s boundaries", data->s2);

    return true;
}
