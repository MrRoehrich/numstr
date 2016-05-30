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
    char buffer [33];

    strcpy(savePath, "OutParaview");
    strcat(savePath, OS_SEP);
    strcat(savePath, "dragonfly");
    sprintf(buffer,"%i",curIter);
    strcat(savePath, buffer);
    strcat(savePath, ".vtk");
    if(!saveDataVtk(data, savePath, curIter))
        return false;

    return true;
}

//--------------------vtk/paraview output----------------------
bool saveDataVtk(const sData* data, const char* vtkFilePath, int curIter)
{
    char buffer [33];

    std::ofstream resultFile(vtkFilePath);
    if(!resultFile) {
        return false;
    }
    resultFile.clear();

    resultFile << "# vtk DataFile Version 3.0" << std::endl;
    resultFile << "vtk output" << curIter << std::endl;
    resultFile << "ASCII" << std::endl;
    resultFile << "DATASET STRUCTURED_GRID" << std::endl;
    resultFile << "DIMENSIONS " << data->nPointsX << " " << data->nPointsY << " 1" << std::endl;
    resultFile << "POINTS " << data->nPoints << " float" << std::endl;
    for(int pointId = 0; pointId < data->nPoints; pointId++) {
        resultFile << data->points[pointId].x << " " << data->points[pointId].y << " "
                   << "0" << std::endl;
    }
    resultFile << std::endl;

    resultFile << "CELL_DATA " << data->nCells << std::endl;
    resultFile << "SCALARS "
               << "phi "
               << "float" << std::endl;
    resultFile << "LOOKUP_TABLE dafault " << std::endl;
    for(int i = 0; i < data->nCells; i++) {
        resultFile << data->cells[i].phi << std::endl;
    }
    resultFile << std::endl;
    resultFile.close();
    return true;
}
