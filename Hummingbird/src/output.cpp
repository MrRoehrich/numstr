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

#ifdef WIN32
#define OS_SEP "\\"
#else
#define OS_SEP "/"
#endif

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "data.h"
#include "output.h"

//-----------------------------------------------------

bool output(sData* data)
{
    char savePath[80];

    if(!computeVelocity(data)) {
        return false;
    }

    if(!computePressureCoef(data)) {
        return false;
    }

    std::cout << "\nOutput:\n-------\n";

    showScalar(data, "scalar 1");

    strcpy(savePath, "OutParaview");
    strcat(savePath, OS_SEP);
    strcat(savePath, "hummingbird.vtk");
    if(!saveVtk(data, savePath, data->s1)) {
        return false;
    }

    strcpy(savePath, "OutMatlab");
    strcat(savePath, OS_SEP);
    strcat(savePath, "grid");
    if(!saveMatlabGrid(data, savePath)) {
        return false;
    }

    strcpy(savePath, "OutMatlab");
    strcat(savePath, OS_SEP);
    strcat(savePath, "scalar.dat");
    if(!saveMatlabScalar(data, savePath, data->s1)) {
        return false;
    }

    return true;
}

//--- compute the velocity field -------
bool computeVelocity(sData* data)
{
    // compute velocity field for scalar s1
    for(int i = 1; i < data->nX - 1; i++) {
        for(int j = 1; j < data->nY - 1; j++) {
            // calculate derivatives at inner grid points
            data->ds1dxi[i][j] = (data->s1[i + 1][j] - data->s1[i - 1][j]) / 2.;
            data->ds1det[i][j] = (data->s1[i][j + 1] - data->s1[i][j - 1]) / 2.;
            data->u[i][j] = data->ds1dxi[i][j] * data->dxidx[i][j] + data->ds1det[i][j] * data->detdx[i][j];
            data->v[i][j] = data->ds1dxi[i][j] * data->dxidy[i][j] + data->ds1det[i][j] * data->detdy[i][j];
        }
        // at lower boundary
        data->ds1dxi[i][0] = (data->s1[i + 1][0] - data->s1[i - 1][0]) / 2.;
        data->ds1det[i][0] = (data->s1[i][1] - data->s1[i][0]);
        data->u[i][0] = data->ds1dxi[i][0] * data->dxidx[i][0] + data->ds1det[i][0] * data->detdx[i][0];
        // at  upper boundary
        data->ds1dxi[i][data->nY - 1] = (data->s1[i + 1][data->nY - 1] - data->s1[i - 1][data->nY - 1]) / 2.;
        data->ds1det[i][data->nY - 1] = (data->s1[i][data->nY - 1] - data->s1[i][data->nY - 2]);
        data->u[i][data->nY - 1] = data->ds1dxi[i][data->nY - 1] * data->dxidx[i][data->nY - 1] +
            data->ds1det[i][data->nY - 1] * data->detdx[i][data->nY - 1];

        // at lower boundary
        data->ds1dxi[i][0] = (data->s1[i + 1][0] - data->s1[i - 1][0]) / 2.;
        data->ds1det[i][0] = (data->s1[i][1] - data->s1[i][0]);
        data->v[i][0] = data->ds1dxi[i][0] * data->dxidy[i][0] + data->ds1det[i][0] * data->detdy[i][0];
        // at  upper boundary
        data->ds1dxi[i][data->nY - 1] = (data->s1[i + 1][data->nY - 1] - data->s1[i - 1][data->nY - 1]) / 2.;
        data->ds1det[i][data->nY - 1] = (data->s1[i][data->nY - 1] - data->s1[i][data->nY - 2]);
        data->v[i][data->nY - 1] = data->ds1dxi[i][data->nY - 1] * data->dxidy[i][data->nY - 1] +
            data->ds1det[i][data->nY - 1] * data->detdy[i][data->nY - 1];
    }
    for(int j = 1; j < data->nY - 1; ++j) {
        // at left boundary
        data->ds1dxi[0][j] = (data->s1[1][j] - data->s1[0][j]);
        data->ds1det[0][j] = (data->s1[0][j + 1] - data->s1[0][j - 1]) / 2.;
        data->v[0][j] = data->ds1dxi[0][j] * data->dxidy[0][j] + data->ds1det[0][j] * data->detdy[0][j];
        // at  right boundary
        data->ds1dxi[data->nX - 1][j] = (data->s1[data->nX - 1][j] - data->s1[data->nX - 2][j]);
        data->ds1det[data->nX - 1][j] = (data->s1[data->nX - 1][j + 1] - data->s1[data->nX - 1][j - 1]) / 2.;
        data->v[data->nX - 1][j] = data->ds1dxi[data->nX - 1][j] * data->dxidy[data->nX - 1][j] +
            data->ds1det[data->nX - 1][j] * data->detdy[data->nX - 1][j];

        // at left boundary
        data->ds1dxi[0][j] = (data->s1[1][j] - data->s1[0][j]);
        data->ds1det[0][j] = (data->s1[0][j + 1] - data->s1[0][j - 1]) / 2.;
        data->u[0][j] = data->ds1dxi[0][j] * data->dxidx[0][j] + data->ds1det[0][j] * data->detdx[0][j];
        // at  right boundary
        data->ds1dxi[data->nX - 1][j] = (data->s1[data->nX - 1][j] - data->s1[data->nX - 2][j]);
        data->ds1det[data->nX - 1][j] = (data->s1[data->nX - 1][j + 1] - data->s1[data->nX - 1][j - 1]) / 2.;
        data->u[data->nX - 1][j] = data->ds1dxi[data->nX - 1][j] * data->dxidx[data->nX - 1][j] +
            data->ds1det[data->nX - 1][j] * data->detdx[data->nX - 1][j];
    }

    data->ds1dxi[0][0] = (data->s1[1][0] - data->s1[0][0]);
    data->ds1det[0][0] = (data->s1[0][1] - data->s1[0][0]);
    data->u[0][0] = data->ds1dxi[0][0] * data->dxidx[0][0] + data->ds1det[0][0] * data->detdx[0][0];
    data->v[0][0] = data->ds1dxi[0][0] * data->dxidy[0][0] + data->ds1det[0][0] * data->detdy[0][0];

    data->ds1dxi[data->nX - 1][0] = (data->s1[data->nX - 1][0] - data->s1[data->nX - 2][0]);
    data->ds1det[data->nX - 1][0] = (data->s1[data->nX - 1][1] - data->s1[data->nX - 1][0]);
    data->u[data->nX - 1][0] = data->ds1dxi[data->nX - 1][0] * data->dxidx[data->nX - 1][0] +
        data->ds1det[data->nX - 1][0] * data->detdx[data->nX - 1][0];
    data->v[data->nX - 1][0] = data->ds1dxi[data->nX - 1][0] * data->dxidy[data->nX - 1][0] +
        data->ds1det[data->nX - 1][0] * data->detdy[data->nX - 1][0];

    data->ds1dxi[0][data->nY - 1] = (data->s1[1][data->nY - 1] - data->s1[0][data->nY - 1]);
    data->ds1det[0][data->nY - 1] = (data->s1[0][data->nY - 1] - data->s1[0][data->nY - 2]);
    data->u[0][data->nY - 1] = data->ds1dxi[0][data->nY - 1] * data->dxidx[0][data->nY - 1] +
        data->ds1det[0][data->nY - 1] * data->detdx[0][data->nY - 1];
    data->v[0][data->nY - 1] = data->ds1dxi[0][data->nY - 1] * data->dxidy[0][data->nY - 1] +
        data->ds1det[0][data->nY - 1] * data->detdy[0][data->nY - 1];

    data->ds1dxi[data->nX - 1][data->nY - 1] =
        (data->s1[data->nX - 1][data->nY - 1] - data->s1[data->nX - 2][data->nY - 1]);
    data->ds1det[data->nX - 1][data->nY - 1] =
        (data->s1[data->nX - 1][data->nY - 1] - data->s1[data->nX - 1][data->nY - 2]);

    data->u[data->nX - 1][data->nY - 1] =
        data->ds1dxi[data->nX - 1][data->nY - 1] * data->dxidx[data->nX - 1][data->nY - 1] +
        data->ds1det[data->nX - 1][data->nY - 1] * data->detdx[data->nX - 1][data->nY - 1];
    data->v[data->nX - 1][data->nY - 1] =
        data->ds1dxi[data->nX - 1][data->nY - 1] * data->dxidy[data->nX - 1][data->nY - 1] +
        data->ds1det[data->nX - 1][data->nY - 1] * data->detdy[data->nX - 1][data->nY - 1];

    return true;
}

//--- compute the pressure coefficient -------
bool computePressureCoef(sData* data)
{
    // compute pressure coefficent for scalar s1
    double vRefSquared;
    if(data->potentialFunc == CHANNELFLOW) {
        vRefSquared = data->uIn * data->uIn;
    } else if(data->potentialFunc == PARALLELFLOW) {
        vRefSquared = (data->uInfty * data->uInfty) + (data->vInfty * data->vInfty);
    } else if(data->potentialFunc == STAGNATION_POINT) {
        vRefSquared = MAXDOUBLE;
    } else if(data->potentialFunc == SOURCE) {
        vRefSquared = MAXDOUBLE;
    } else if(data->potentialFunc == POTENTIAL_VORTEX) {
        vRefSquared = MAXDOUBLE;
    } else if (data->potentialFunc == DIPOL) {
        vRefSquared = MAXDOUBLE;
    }
    for(int i = 0; i < data->nX; ++i)
        for(int j = 0; j < data->nY; ++j)
            data->cp[i][j] = 1 - (data->u[i][j] * data->u[i][j] + data->v[i][j] * data->v[i][j]) / vRefSquared;
    return true;
}

void printValue(const sData* data, const char* scalarName, double** value)
{
    const int maxHoriz = 5;
    const int maxVert = 5;

    std::cout.precision(2);

    std::cout << "\net\t------------------------- " << scalarName << " -------------------------\n"
              << "^\n"
              << "|\n";

    double iStep, jStep;
    if(data->nX < maxVert) {
        iStep = 1;
    } else {
        iStep = data->nX / (double)maxVert;
    }
    if(data->nY < maxHoriz) {
        jStep = 1;
    } else {
        jStep = data->nY / (double)maxHoriz;
    }

    double i, j = data->nY - 1 + jStep;
    while(j > 0) {
        j -= jStep;
        if(j < 1) {
            j = 0;
        }
        std::cout << std::fixed << (int)j << "\t";

        i = -iStep;
        while(i < data->nX - 1) {
            i += iStep;
            if(i > data->nX - 2) {
                i = data->nX - 1;
            }
            std::cout.setf(std::ios::showpos);
            std::cout << std::scientific << value[(int)i][(int)j] << "  ";
            std::cout.unsetf(std::ios::showpos);
        }
        std::cout << "\n|\n";
    }
    std::cout << " --\t";

    i = -iStep;
    while(i < data->nX - 1) {
        i += iStep;
        if(i > data->nX - 2)
            i = data->nX - 1;
        std::cout << "   -" << (int)i << "-    ";
    }
    std::cout << " ->xi\n\n";
}

//--- shows the scalar on the console------------------
void showScalar(const sData* data, const char* scalarName)
{
    printValue(data, "scalar 1", data->s1);
    printValue(data, "scalar 2", data->s2);
    printValue(data, "u", data->u);
    printValue(data, "v", data->v);
    printValue(data, "cp", data->cp);
}

//--- save the scalar on Matlab readable format -------
bool saveMatlabScalar(const sData* data, const char* scalarFilePath, double** s)
{
    std::ofstream scalarFile(scalarFilePath);
    if(!scalarFile) {
        return false;
    }
    scalarFile.clear();
    for(int i = 0; i < data->nX; i++) {
        for(int j = 0; j < data->nY; j++) {
            scalarFile << " " << s[i][j];
        }
        scalarFile << std::endl;
    }
    scalarFile.close();

    return true;
}

//--- saves the grid in Matlab readable format --------
bool saveMatlabGrid(const sData* data, const char* gridFilePath)
{
    char fileNameX[80];
    char fileNameY[80];

    sprintf(fileNameX, "%sX.dat", gridFilePath);
    sprintf(fileNameY, "%sY.dat", gridFilePath);
    std::ofstream meshXFile(fileNameX);
    std::ofstream meshYFile(fileNameY);
    if(!meshXFile) {
        return false;
    }
    if(!meshYFile) {
        return false;
    }
    meshXFile.clear();
    meshYFile.clear();
    for(int i = 0; i < data->nX; i++) {
        for(int j = 0; j < data->nY; j++) {
            meshXFile << data->x[i][j] << " ";
            meshYFile << data->y[i][j] << " ";
        }
        meshXFile << std::endl;
        meshYFile << std::endl;
    }
    meshXFile.close();
    meshYFile.close();

    return true;
}

//--- save grid and scalar in vtk format --------------
bool saveVtk(const sData* data, const char* vtkFilePath, double** s)
{

    // open file
    std::ofstream dataFile(vtkFilePath);
    if(!dataFile) {
        return false;
    }
    dataFile.clear();

    // write header
    dataFile << "# vtk DataFile Version 2.0" << std::endl;
    dataFile << "Results formatted in VTK File Format" << std::endl;
    dataFile << "ASCII" << std::endl;
    dataFile << "DATASET STRUCTURED_GRID" << std::endl;
    dataFile << "DIMENSIONS " << data->nX << " " << data->nY << " 1" << std::endl;
    dataFile << "POINTS " << (data->nX) * (data->nY) << " float" << std::endl;

    // write point coordinates
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->x[i][j] << " " << data->y[i][j] << " 0" << std::endl;
        }
    }
    dataFile << std::endl;

    // write the scalar data to result file
    dataFile << "POINT_DATA " << (data->nX) * (data->nY) << std::endl;
    dataFile << "SCALARS scalar1 float" << std::endl;
    dataFile << "LOOKUP_TABLE default" << std::endl;
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->s1[i][j] << std::endl;
        }
    }
    dataFile << std::endl;

    // write the scalar data to result file
    dataFile << "SCALARS scalar2 float" << std::endl;
    dataFile << "LOOKUP_TABLE default" << std::endl;
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->s2[i][j] << std::endl;
        }
    }
    dataFile << std::endl;

    // write the scalar data to result file
    dataFile << "SCALARS  u float" << std::endl;
    dataFile << "LOOKUP_TABLE default" << std::endl;
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->u[i][j] << std::endl;
        }
    }
    dataFile << std::endl;

    dataFile << "SCALARS  v float" << std::endl;
    dataFile << "LOOKUP_TABLE default" << std::endl;
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->v[i][j] << std::endl;
        }
    }
    dataFile << std::endl;

    dataFile << "SCALARS  cp float" << std::endl;
    dataFile << "LOOKUP_TABLE default" << std::endl;
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->cp[i][j] << std::endl;
        }
    }
    dataFile << std::endl;

    dataFile << "VECTORS  velocity float" << std::endl;
    for(int j = 0; j < data->nY; j++) {
        for(int i = 0; i < data->nX; i++) {
            dataFile << data->u[i][j] << " ";
            dataFile << data->v[i][j] << " ";
            dataFile << 0. << std::endl;
        }
    }

    dataFile.close();

    return true;
}
