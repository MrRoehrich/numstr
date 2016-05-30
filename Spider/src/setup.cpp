/***************************************************************************
 *   Copyright (C) 2006-2014 by  Institute of Combustion Technology        *
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
    std::cout << "\nSetup:\n-------\n";

    data->nPointsX = data->nCellsX + 1;
    data->nPointsY = data->nCellsY + 1;

    // allocate memory
    data->pointsX = new double*[data->nPointsX];
    data->pointsY = new double*[data->nPointsX];
    for(int i = 0; i < data->nPointsX; i++) {
        data->pointsX[i] = new double[data->nPointsY];
        data->pointsY[i] = new double[data->nPointsY];
    }

    // calculate point coordinates
    double deltaX = (data->xMax - data->xMin) / data->nCellsX;
    double deltaY = (data->yMax - data->yMin) / data->nCellsY;
    for(int i = 0; i < data->nPointsX; i++) {
        for(int j = 0; j < data->nPointsY; j++) {
            data->pointsX[i][j] = data->xMin + deltaX * i;
            data->pointsY[i][j] = data->yMin + deltaY * j;
        }
    }

    return true;
}
