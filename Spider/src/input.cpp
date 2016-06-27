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
#include <string>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "data.h"
#include "input.h"

//------------------------------------------------------
bool input(const char* cfgFilePath, sData* data)
{
    char line[256] = " ", token[16] = " ";
    std::cout << "\nInput:\n-------\n";

    //////////////////////
    // READ CFG FILE    //
    //////////////////////

    // open input file
    std::ifstream cfgFile(cfgFilePath);
    if(!cfgFile) {
        return false;
    }

    // read input file line by line
    while(!cfgFile.eof()) {
        cfgFile.getline(line, 255);
        if(sscanf(line, "%15s", token) <= 0) {
            continue;
        };

        if(!strcmp(token, "#")) {
            // skip comment lines
        } else if(!strcmp(token, "nCellsX")) {
            sscanf(line, "%15s %d", token, &data->nCellsX);
        } else if(!strcmp(token, "nCellsY")) {
            sscanf(line, "%15s %d", token, &data->nCellsY);
        } else if(!strcmp(token, "xMin")) {
            sscanf(line, "%15s %lf", token, &data->xMin);
        } else if(!strcmp(token, "xMax")) {
            sscanf(line, "%15s %lf", token, &data->xMax);
        } else if(!strcmp(token, "yMin")) {
            sscanf(line, "%15s %lf", token, &data->yMin);
        } else if(!strcmp(token, "yMax")) {
            sscanf(line, "%15s %lf", token, &data->yMax);
        } else if(!strcmp(token, "scalarTypeN")) {
            sscanf(line, "%15s %d", token, &data->scalarTypeN);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token)) {
                if(data->scalarTypeN == 2) {
                    if(!strcmp(token, "valueN"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueNX, &data->valueNY);
                } else if(data->scalarTypeN == 1)
                    sscanf(line, "%15s %lf", token, &data->valueN);
            }
        } else if(!strcmp(token, "scalarTypeE")) {
            sscanf(line, "%15s %d", token, &data->scalarTypeE);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token)) {
                if(data->scalarTypeE == 2) {
                    if(!strcmp(token, "valueE"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueEX, &data->valueEY);
                } else if(data->scalarTypeE == 1)
                    sscanf(line, "%15s %lf", token, &data->valueE);
            }
        } else if(!strcmp(token, "scalarTypeS")) {
            sscanf(line, "%15s %d", token, &data->scalarTypeS);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token)) {
                if(data->scalarTypeS == 2) {
                    if(!strcmp(token, "valueS"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueSX, &data->valueSY);
                } else if(data->scalarTypeS == 1)
                    sscanf(line, "%15s %lf", token, &data->valueS);
            }
        } else if(!strcmp(token, "scalarTypeW")) {
            sscanf(line, "%15s %d", token, &data->scalarTypeW);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token)) {
                if(data->scalarTypeW == 2) {
                    if(!strcmp(token, "valueW"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueWX, &data->valueWY);
                } else if(data->scalarTypeW == 1)
                    sscanf(line, "%15s %lf", token, &data->valueW);
            }
        } else if(!strcmp(token, "velocityTypeN")) {
            sscanf(line, "%15s %d", token, &data->velocityTypeN);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token))
                if(data->velocityTypeN == 1)
                    if(!strcmp(token, "valueN"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueNU, &data->valueNV);
        } else if(!strcmp(token, "velocityTypeE")) {
            sscanf(line, "%15s %d", token, &data->velocityTypeE);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token))
                if(data->velocityTypeE == 1)
                    if(!strcmp(token, "valueE"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueEU, &data->valueEV);
        } else if(!strcmp(token, "velocityTypeS")) {
            sscanf(line, "%15s %d", token, &data->velocityTypeS);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token))
                if(data->velocityTypeS == 1)
                    if(!strcmp(token, "valueS"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueSU, &data->valueSV);
        } else if(!strcmp(token, "velocityTypeW")) {
            sscanf(line, "%15s %d", token, &data->velocityTypeW);
            cfgFile.getline(line, 255);
            if(sscanf(line, "%15s", token))
                if(data->velocityTypeW == 1)
                    if(!strcmp(token, "valueW"))
                        sscanf(line, "%15s %lf %lf", token, &data->valueWU, &data->valueWV);
        } else if(!strcmp(token, "initConds")) {
            sscanf(line, "%15s %lf", token, &data->initialCondition);
        } else {
            std::cout << "unknown token: " << token << std::endl;
        }
    }
    cfgFile.close();

    return true;
}