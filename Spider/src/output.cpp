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

    // scalar cell boundary conditions
    int cId;
    meshFile << std::endl;
    meshFile << "boundaryConditionsScalar" << std::endl;
    meshFile << "# type: 1=DIRICHLET" << std::endl;
    meshFile << "# type: 2=NEUMANN" << std::endl;
    meshFile << "# cellId type value" << std::endl;
    for(int j = 0; j < data->nCellsY; ++j)
        for(int i = 0; i < data->nCellsX; ++i) {
            cId = i + j * data->nCellsX;
            if(cId < data->nCellsX) { // unterer Rand
                if(data->scalarTypeS == 1)
                    meshFile << cId << " " << data->scalarTypeS << " " << data->valueS << std::endl;
                else if(data->scalarTypeS == 2) {
                    if((cId + data->nCellsX) % data->nCellsX == 0 && data->scalarTypeW == 2) { // linker Rand
                        meshFile << cId << " " << data->scalarTypeS << " " << data->valueWX << " " << data->valueSY
                                 << std::endl;
                        continue;
                    } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0 &&
                        data->scalarTypeE == 2) { // rechter Rand
                        meshFile << cId << " " << data->scalarTypeS << " " << data->valueEX << " " << data->valueSY
                                 << std::endl;
                        continue;
                    } else
                        meshFile << cId << " " << data->scalarTypeS << " " << data->valueSX << " " << data->valueSY
                                 << std::endl;
                }
            } else if((data->nCellsX * data->nCellsY) - (cId + 1) < data->nCellsX && data->scalarTypeN) { // oberer Rand
                if(data->scalarTypeN == 1)
                    meshFile << cId << " " << data->scalarTypeN << " " << data->valueN << std::endl;
                else if(data->scalarTypeN == 2) {
                    if((cId + data->nCellsX) % data->nCellsX == 0 && data->scalarTypeW == 2) { // linker Rand
                        meshFile << cId << " " << data->scalarTypeN << " " << data->valueWX << " " << data->valueNY
                                 << std::endl;
                        continue;
                    } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0 &&
                        data->scalarTypeE == 2) { // rechter Rand
                        meshFile << cId << " " << data->scalarTypeN << " " << data->valueEX << " " << data->valueNY
                                 << std::endl;
                        continue;
                    } else
                        meshFile << cId << " " << data->scalarTypeN << " " << data->valueNX << " " << data->valueNY
                                 << std::endl;
                }
            } else if((cId + data->nCellsX) % data->nCellsX == 0) { // linker Rand
                if(data->scalarTypeW == 1) {
                    meshFile << cId << " " << data->scalarTypeW << " " << data->valueW << std::endl;
                } else if(data->scalarTypeW == 2) {
                    meshFile << cId << " " << data->scalarTypeN << " " << data->valueWX << " " << data->valueNY
                             << std::endl;
                    continue;
                }
            } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0) { // rechter Rand
                if(data->scalarTypeE == 1) {
                    meshFile << cId << " " << data->scalarTypeE << " " << data->valueE << std::endl;
                } else if(data->scalarTypeE == 2) {
                    meshFile << cId << " " << data->scalarTypeN << " " << data->valueEX << " " << data->valueNY
                             << std::endl;
                    continue;
                }
            }
        }
    meshFile << std::endl;

    // cell boundary conditions
    meshFile << std::endl;
    meshFile << "boundaryConditionsVelocity" << std::endl;
    meshFile << "# type: 1=DIRICHLET" << std::endl;
    meshFile << "# type: 2=NEUMANN" << std::endl;
    meshFile << "# cellId type value" << std::endl;
    for(int j = 0; j < data->nCellsY; ++j)
        for(int i = 0; i < data->nCellsX; ++i) {
            cId = i + j * data->nCellsX;
            if(cId < data->nCellsX) { // unterer Rand
                if(data->velocityTypeS == 1) {
                    if((cId + data->nCellsX) % data->nCellsX == 0 && data->velocityTypeW == 1) { // linker Rand
                        meshFile << cId << " " << data->velocityTypeS << " " << data->valueWU << " " << data->valueSV
                                 << std::endl;
                        continue;
                    } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0 &&
                        data->velocityTypeE == 1) { // rechter Rand
                        meshFile << cId << " " << data->velocityTypeS << " " << data->valueEU << " " << data->valueSV
                                 << std::endl;
                        continue;
                    } else
                        meshFile << cId << " " << data->velocityTypeS << " " << data->valueSU << " " << data->valueSV
                                 << std::endl;
                }
            } else if((data->nCellsX * data->nCellsY) - (cId + 1) < data->nCellsX &&
                data->velocityTypeN) { // oberer Rand
                if(data->velocityTypeN == 1) {
                    if((cId + data->nCellsX) % data->nCellsX == 0 && data->velocityTypeN == 1) { // linker Rand
                        meshFile << cId << " " << data->velocityTypeN << " " << data->valueWU << " " << data->valueNV
                                 << std::endl;
                        continue;
                    } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0 &&
                        data->velocityTypeE == 1) { // rechter Rand
                        meshFile << cId << " " << data->velocityTypeN << " " << data->valueEU << " " << data->valueNV
                                 << std::endl;
                        continue;
                    } else
                        meshFile << cId << " " << data->velocityTypeN << " " << data->valueNU << " " << data->valueNV
                                 << std::endl;
                }
            } else if((cId + data->nCellsX) % data->nCellsX == 0 && data->velocityTypeW == 1) { // linker Rand
                meshFile << cId << " " << data->velocityTypeN << " " << data->valueWU << " " << data->valueWV
                         << std::endl;
                continue;
            } else if((cId - (data->nCellsX - 1)) % data->nCellsX == 0 && data->velocityTypeE == 1) { // rechter Rand
                meshFile << cId << " " << data->velocityTypeN << " " << data->valueEU << " " << data->valueEV
                         << std::endl;
                continue;
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
    meshFile << std::endl;

    // sources
    meshFile << "sources" << std::endl;
    meshFile << "# cellId value" << std::endl;

    meshFile.close();

    return true;
}
