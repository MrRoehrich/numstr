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
#include <fstream>
#include <iostream>
#include <stdio.h>

#include "data.h"
#include "output.h"

//------------------------------------------------------
bool output(sData* data)
{
    std::cout << "\nOutput:\n-------\n";

    if(!saveMesh(data)) {
        return false;
    }

    return true;
}

//------------------------------------------------------
bool saveMesh(const sData* data)
{
    const char* fileName = "dragonfly.mesh";
    int id = 0;

    std::ofstream meshFile(fileName);

    if(!meshFile) {
        return false;
    }
    meshFile.clear();

    // point dimensions
    meshFile << "pointDimensions" << std::endl;
    meshFile << data->nPointsX << " " << data->nPointsY << std::endl;

    // point coordinates
    meshFile << std::endl;
    meshFile << "points" << std::endl;
    meshFile << "# id x y" << std::endl;
    for(int j = 0; j < data->nPointsY; j++) {
        for(int i = 0; i < data->nPointsX; i++) {
            meshFile << id << " " << data->pointsX[i][j] << " " << data->pointsY[i][j] << std::endl;
            id++;
        }
    }

    // cell boundary conditions
    int cId;
    meshFile << std::endl;
    meshFile << "boundaryConditions" << std::endl;
    meshFile << "# type: 1=DIRICHLET" << std::endl;
    meshFile << "# type: 2=NEUMANN" << std::endl;
    meshFile << "# cellId type value" << std::endl;
    for(int j = 0; j < data->nCellsY; ++j)
        for(int i = 0; i < data->nCellsX; ++i) {
            cId = i + j * data->nCellsX;
            if(cId < data->nCellsX && data->typeS) {
                meshFile << cId << " " << data->typeS << " " << data->valueS << std::endl;
                continue;
            } else if((data->nCellsX * data->nCellsY) - (cId + 1) < data->nCellsX && data->typeN) {
                meshFile << cId << " " << data->typeN << " " << data->valueN << std::endl;
                continue;
            }
            if((cId + data->nCellsX) % data->nCellsX == 0 && data->typeW) {
                meshFile << cId << " " << data->typeW << " " << data->valueW << std::endl;
            } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0 && data->typeE) {
                meshFile << cId << " " << data->typeE << " " << data->valueE << std::endl;
            }
        }
        meshFile << std::endl;

    // initial conditions conditions
    meshFile << "initialConditions" << std::endl;
    meshFile << "# cellId value" << std::endl;
    for(int j = 0; j < data->nCellsY; ++j)
        for(int i = 0; i < data->nCellsX; ++i) {
            cId = i + j * data->nCellsX;
            meshFile << cId << " " << data->initialCondition << std::endl;
        }

    meshFile.close();

    return true;
}
