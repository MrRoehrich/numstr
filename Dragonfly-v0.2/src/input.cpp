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
#include <fstream>
#include <iostream>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "data.h"
#include "input.h"

//------------------------------------------------------
bool input(const char* cfgFilePath, const char* meshFilePath, sData* data)
{

    int section, lineNo;
    sCell* curCell = 0;
    sPoint* curPoint = 0;
    char line[256] = " ";
    char token[32] = " ";
    int pointId;
    double bValue;
    double bValueX;
    double bValueY;
    int bType;
    double x;
    double y;
    int cellId;
    double initPhi;
    double source;

    //////////////////////
    // READ JOB FILE  //
    //////////////////////

    // open input file
    std::ifstream cfgFile(cfgFilePath);
    if(!cfgFile) {
        return false;
    }
    // read input file line by line
    lineNo = 0;
    while(!cfgFile.eof()) {
        lineNo++;
        cfgFile.getline(line, 255);
        if(sscanf(line, "%15s", token) < 1) {
            continue;
        };

        if(!strcmp(token, "#")) {
            // skip comment lines
            // numerical settings
        } else if(!strcmp(token, "maxTime")) {
            if(sscanf(line, "%15s %lf", token, &data->maxTime) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
        } else if(!strcmp(token, "solverMethod")) {
            if(sscanf(line, "%15s %d", token, &data->solverMethod) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
        } else if(!strcmp(token, "solverType")) {
            if(sscanf(line, "%15s %d", token, &data->solverType) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
        } else if(!strcmp(token, "maxIter")) {
            if(sscanf(line, "%15s %d", token, &data->maxIter) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
        } else if(!strcmp(token, "residuum")) {
            if(sscanf(line, "%15s %lf", token, &data->residuum) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
            // physical settings
        } else if(!strcmp(token, "alpha")) {
            if(sscanf(line, "%15s %lf", token, &data->alpha) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
        } else if(!strcmp(token, "rho")) {
            if(sscanf(line, "%15s %lf", token, &data->rho) != 2) {
                return error(cfgFilePath, lineNo, line);
            };
        } else if(!strcmp(token, "velocity")) {
            if(sscanf(line, "%15s %lf %lf", token, &data->u, &data->v) != 3) {
                return error(cfgFilePath, lineNo, line);
            };
        } else {
            std::cout << "unknown token: " << token << std::endl;
            return error(cfgFilePath, lineNo, line);
        }
    }
    cfgFile.close();

    /////////////////////
    // READ FACE FILE  //
    /////////////////////

    // open face file
    std::ifstream meshFile(meshFilePath);
    if(!meshFile) {
        return false;
    }

    // read face data
    section = 0;
    lineNo = 0;
    int nX, nY;
    while(!meshFile.eof()) {
        lineNo++;
        meshFile.getline(line, 255);
        if(sscanf(line, "%31s", token) < 1) {
            continue;
        };

        if(!strcmp(token, "#")) {
            // skip comment lines
        } else if(!strcmp(token, "pointDimensions")) {
            if(sscanf(line, "%31s", token) != 1) {
                return error(meshFilePath, lineNo, line);
            };
            section = 1;
        } else if(!strcmp(token, "points")) {
            if(sscanf(line, "%31s", token) != 1) {
                return error(meshFilePath, lineNo, line);
            };
            section = 2;
        } else if(!strcmp(token, "boundaryConditionsScalar")) {
            if(sscanf(line, "%31s", token) != 1) {
                return error(meshFilePath, lineNo, line);
            };
            section = 3;
        } else if(!strcmp(token, "boundaryConditionsVelocity")) {
            if(sscanf(line, "%31s", token) != 1) {
                return error(meshFilePath, lineNo, line);
            };
            section = 4;
        } else if(!strcmp(token, "initialConditions")) {
            if(sscanf(line, "%31s", token) != 1) {
                return error(meshFilePath, lineNo, line);
            };
            section = 5;
        } else if(!strcmp(token, "sources")) {
            if(sscanf(line, "%31s", token) != 1) {
                return error(meshFilePath, lineNo, line);
            };
            section = 6;
        } else if(section == 1) { // reading dimensions
            if(sscanf(line, "%d %d", &nX, &nY) != 2) {
                return error(meshFilePath, lineNo, line);
            };
            data->nPointsX = nX;
            data->nPointsY = nY;
            data->nCellsX = data->nPointsX - 1;
            data->nCellsY = data->nPointsY - 1;
            data->nPoints = nX * nY;
            data->nCells = (nX - 1) * (nY - 1);
            data->nFaces = 2 * nX * nY - nX - nY;
            data->points = new sPoint[data->nPoints];
            data->cells = new sCell[data->nCells];
            for(int cId = 0; cId < data->nCells; cId++) {
                data->cells[cId].id = cId;
                data->cells[cId].p = 0. + (cId - (data->nCellsX - 1)) % (data->nCellsX);
            }
            data->faces = new sFace[data->nFaces];
        } else if(section == 2) { // reading point data section
            if(sscanf(line, "%d %lf %lf", &pointId, &x, &y) != 3) {
                return error(meshFilePath, lineNo, line);
            }
            curPoint = &data->points[pointId];
            curPoint->id = pointId;
            curPoint->x = x;
            curPoint->y = y;
        } else if(section == 3) { // reading boundary condition section
            if(sscanf(line, "%d %d %lf %lf", &cellId, &bType, &bValueX, &bValueY) == 4) {
                curCell = &data->cells[cellId];
                // curCell->id = cellId;
                curCell->bTypeScalar = bType;
                curCell->bValueScalarX = bValueX;
                curCell->bValueScalarY = bValueY;

            } else if(sscanf(line, "%d %d %lf", &cellId, &bType, &bValue) != 3)
                return error(meshFilePath, lineNo, line);
            else {
                curCell = &data->cells[cellId];
                curCell->bTypeScalar = bType;
                curCell->bValueScalar = bValue;
            }
        } else if(section == 4) { // reading boundary condition section
            if(sscanf(line, "%d %d %lf %lf", &cellId, &bType, &bValueX, &bValueY) == 4) {
                curCell = &data->cells[cellId];
                curCell->bTypeVelocity = bType;
                curCell->bValueU = bValueX;
                curCell->bValueV = bValueY;

            } else if(sscanf(line, "%d %d", &cellId, &bType) != 2) {
                return error(meshFilePath, lineNo, line);
            } else {
                curCell = &data->cells[cellId];
                curCell->bTypeVelocity = bType;
            }
        } else if(section == 5) { // reading initial condition section
            if(sscanf(line, "%d %lf", &cellId, &initPhi) != 2) {
                return error(meshFilePath, lineNo, line);
            }
            curCell = &data->cells[cellId];
            curCell->phi = initPhi;
        } else if(section == 6) { // reading source terms
            if(sscanf(line, "%d %lf", &cellId, &source) != 2) {
                return error(meshFilePath, lineNo, line);
            } else {
                curCell = &data->cells[cellId];
                curCell->s = source;
            }
        }
    }
    meshFile.close();
    return true;
}

//------------------------------------------------------
bool error(const char* filePath, int lineNo, const char* line)
{
    std::cout << "ERROR reading " << filePath << ", line " << lineNo << ":" << std::endl;
    std::cout << "\t" << line << std::endl << std::endl;
    return false;
}
