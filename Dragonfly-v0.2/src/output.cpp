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

#ifdef WIN32
#define OS_SEP "\\"
#else
#define OS_SEP "/"
#endif

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#include "data.h"
#include "input.h"
#include "output.h"

//------------------------------------------------------
bool output(sData* data, int curIter)
{
    char savePath[80];

    strcpy(savePath, "OutParaview");
    strcat(savePath, OS_SEP);
    strcat(savePath, "dragonfly");
    if(!saveDataVtk(data, savePath, curIter))
        return false;

    return true;
}

//--------------------vtk/paraview output----------------------
bool saveDataVtk(const sData* data, const char* vtkFilePath, int curIter)
{
    char buffer[33];
    char savePath[80];
    strcpy(savePath, vtkFilePath);
    strcat(savePath, "Phi");
    sprintf(buffer, "%i", curIter);
    strcat(savePath, buffer);
    strcat(savePath, ".vtk");
    std::ofstream resultFilePhi(savePath);
    if(!resultFilePhi) {
        return false;
    }
    resultFilePhi.clear();

    resultFilePhi << "# vtk DataFile Version 3.0" << std::endl;
    resultFilePhi << "vtk output" << curIter << std::endl;
    resultFilePhi << "ASCII" << std::endl;
    resultFilePhi << "DATASET STRUCTURED_GRID" << std::endl;
    resultFilePhi << "DIMENSIONS " << data->nPointsX << " " << data->nPointsY << " 1" << std::endl;
    resultFilePhi << "POINTS " << data->nPoints << " float" << std::endl;
    for(int pointId = 0; pointId < data->nPoints; pointId++) {
        resultFilePhi << data->points[pointId].x << " " << data->points[pointId].y << " "
                      << "0" << std::endl;
    }
    resultFilePhi << std::endl;

    resultFilePhi << "CELL_DATA " << data->nCells << std::endl;
    resultFilePhi << "SCALARS "
                  << "phi "
                  << "float" << std::endl;
    resultFilePhi << "LOOKUP_TABLE default " << std::endl;
    for(int i = 0; i < data->nCells; i++) {
        resultFilePhi << data->cells[i].phi << std::endl;
    }
    resultFilePhi << std::endl;

    resultFilePhi << "SCALARS "
                  << "p "
                  << "float" << std::endl;
    resultFilePhi << "LOOKUP_TABLE default " << std::endl;
    for(int i = 0; i < data->nCells; i++) {
        resultFilePhi << data->cells[i].p << std::endl;
    }
    resultFilePhi << std::endl;
    
    sCell* curCell;
    double u, v;
    
    resultFilePhi << "SCALARS  u float" << std::endl;
    resultFilePhi << "LOOKUP_TABLE default " << std::endl;
    for(int cId = 0; cId < data->nCells; ++cId) {
        curCell = &data->cells[cId];
        u = (curCell->faces[YP]->u + curCell->faces[XP]->u + curCell->faces[YM]->u + curCell->faces[XM]->u)/4.;
        resultFilePhi << u  << std::endl;
    }
    resultFilePhi << std::endl;
    
    resultFilePhi << "SCALARS  v float" << std::endl;
    resultFilePhi << "LOOKUP_TABLE default " << std::endl;
    for(int cId = 0; cId < data->nCells; ++cId) {
        curCell = &data->cells[cId];
        v = (curCell->faces[YP]->v + curCell->faces[XP]->v + curCell->faces[YM]->v + curCell->faces[XM]->v)/4.;
        resultFilePhi << v  << std::endl;
    }
    resultFilePhi << std::endl;
    
    resultFilePhi << "VECTORS  velocity float" << std::endl;
    for(int cId = 0; cId < data->nCells; ++cId) {
        curCell = &data->cells[cId];
        u = (curCell->faces[YP]->u + curCell->faces[XP]->u + curCell->faces[YM]->u + curCell->faces[XM]->u)/4.;
        v = (curCell->faces[YP]->v + curCell->faces[XP]->v + curCell->faces[YM]->v + curCell->faces[XM]->v)/4.;
        resultFilePhi << u << " ";
        resultFilePhi << v << " ";
        resultFilePhi << 0. << std::endl;
    }
    resultFilePhi << std::endl;

    resultFilePhi.close();

    /*// u
    strcpy(savePath, vtkFilePath);
    strcat(savePath, "U");
    sprintf(buffer, "%i", curIter);
    strcat(savePath, buffer);
    strcat(savePath, ".vtk");
    std::ofstream resultFileU(savePath);
    if(!resultFileU) {
        return false;
    }
    resultFileU.clear();

    resultFileU << "# vtk resultFilePhi Version 3.0" << std::endl;
    resultFileU << "vtk output" << curIter << std::endl;
    resultFileU << "ASCII" << std::endl;
    resultFileU << "DATASET STRUCTURED_GRID" << std::endl;
    resultFileU << "DIMENSIONS " << data->nPointsX << " " << data->nPointsY << " 1" << std::endl;
    resultFileU << "POINTS " << data->nPoints << " float" << std::endl;
    for(int pointId = 0; pointId < data->nPoints; pointId++) {
        resultFileU << data->points[pointId].x << " " << data->points[pointId].y << " "
                    << "0" << std::endl;
    }
    resultFileU << std::endl;

    resultFileU << "CELL_DATA " << data->nCells << std::endl;
    resultFileU << "SCALARS "
                << "u "
                << "float" << std::endl;
    resultFileU << "LOOKUP_TABLE default " << std::endl;
    for(int i = 0; i < data->nCells; i++) {
        resultFileU << data->cells[i].faces[XM]->u << std::endl;
    }
    resultFileU << std::endl;
    resultFileU.close();

    strcpy(savePath, vtkFilePath);
    strcat(savePath, "V");
    sprintf(buffer, "%i", curIter);
    strcat(savePath, buffer);
    strcat(savePath, ".vtk");
    std::ofstream resultFileV(savePath);
    if(!resultFileV) {
        return false;
    }
    resultFileV.clear();

    resultFileV << "# vtk resultFilePhi Version 3.0" << std::endl;
    resultFileV << "vtk output" << curIter << std::endl;
    resultFileV << "ASCII" << std::endl;
    resultFileV << "DATASET STRUCTURED_GRID" << std::endl;
    resultFileV << "DIMENSIONS " << data->nPointsX << " " << data->nPointsY << " 1" << std::endl;
    resultFileV << "POINTS " << data->nPoints << " float" << std::endl;
    for(int pointId = 0; pointId < data->nPoints; pointId++) {
        resultFileV << data->points[pointId].x << " " << data->points[pointId].y << " "
                    << "0" << std::endl;
    }
    resultFileV << std::endl;

    resultFileV << "CELL_DATA " << data->nCells << std::endl;
    resultFileV << "SCALARS "
                << "v "
                << "float" << std::endl;
    resultFileV << "LOOKUP_TABLE default " << std::endl;
    for(int i = 0; i < data->nCells; i++) {
        resultFileV << data->cells[i].faces[YM]->v << std::endl;
    }
    resultFileV << std::endl;
    resultFileV.close();*/
    return true;
}
