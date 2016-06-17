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
#ifndef DATA_H
#define DATA_H

#include <limits>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ABS(x) (((x) > 0) ? (x) : -(x))

#define MAXDOUBLE (std::numeric_limits<double>::max())
#define MINDOUBLE (std::numeric_limits<double>::min())
#define MAXFLOAT (std::numeric_limits<float>::max())
#define MINFLOAT (std::numeric_limits<float>::min())
#define MAXINT (std::numeric_limits<int>::max())
#define MININT (std::numeric_limits<int>::min())

#define INNERCELL 0
#define DIRICHLET 1
#define NEUMANN 2

struct sData;

struct sData {
    int nCellsX;
    int nCellsY;

    int nPointsX;
    int nPointsY;

    double xMin;
    double xMax;
    double yMin;
    double yMax;

    int typeN;
    double valueN;
    double valueNX;
    double valueNY;
    int typeE;
    double valueE;
    double valueEX;
    double valueEY;
    int typeS;
    double valueS;
    double valueSX;
    double valueSY;
    int typeW;
    double valueW;
    double valueWX;
    double valueWY;

    double initialCondition;

    double** pointsX;
    double** pointsY;

    ~sData()
    {
        for(int i = 0; i <= nCellsX; i++) {
            delete pointsX[i];
            delete pointsY[i];
        }
        delete pointsX;
        delete pointsY;
    }
};

#endif
