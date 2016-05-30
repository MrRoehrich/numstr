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
#include <stdio.h>

#include "data.h"
#include "input.h"
#include "output.h"
#include "setup.h"
#include "solve.h"

//------------------------------------------------------
int main(int, char**)
{
    sData* data = new sData;
    const char* cfgFilePath = "dragonfly.in";
    const char* meshFilePath = "dragonfly.mesh";
    // const char* cellFilePath  = "mesh.cells";

    std::cout << "\t Numerische Strömungsmechanik\n"
              << "\t====== FlowSim ======\n";

    // read config from input file
    std::cout << "Input ...\n";
    if(!input(cfgFilePath, meshFilePath, data)) {
        std::cout << "ERROR while reading *.in files ...exiting";
        getchar();
        return 1;
    }

    std::cout << "Setup...\n";
    // setup data (boundaries, initial data, etc.)
    if(!setup(data)) {
        std::cout << "ERROR while data setup...exiting";
        getchar();
        return 1;
    }

    std::cout << "Solve...\n";
    // iterativ solver
    if(data->solverMethod == CALCFLUX) {
        if(!solveCalcFlux(data)) {
            std::cout << "ERROR while solving...exiting";
            getchar();
            return 1;
        }
    } else if(data->solverMethod == PE) {
        if(!solvePe(data)) {
            std::cout << "ERROR while solving...exiting";
            getchar();
            return 1;
        }
    }

    std::cout << "Success...";
    delete data;
    return 0;
}
