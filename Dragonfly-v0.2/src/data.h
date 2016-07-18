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
#ifndef DATA_H
#define DATA_H

#include <limits>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ABS(x) (((x) > 0) ? (x) : -(x))

#define DECIDE(a, b, c) (((a) == (c)) ? (a) : (b))

#define MAXDOUBLE (std::numeric_limits<double>::max())
#define MINDOUBLE (std::numeric_limits<double>::min())
#define MAXINT (std::numeric_limits<int>::max())
#define MININT (std::numeric_limits<int>::min())

#define EPS 1e-10

#define YM 0
#define XP 1
#define YP 2
#define XM 3

#define XMYM 0
#define XPYM 1
#define XPYP 2
#define XMYP 3

#define M 0
#define P 1

// solver method
#define CALCFLUX 0
#define PE 1
#define SIMPLE 2

// solver types
#define CENTRAL 0
#define UPWIND 1
#define HYBRID 2
#define POWER 3
#define EXPONENTIAL 4

// boundary types
#define INNERCELL 0
#define DIRICHLET 1
#define NEUMANN 2

// cell types
struct sFace;
struct sCell;
struct sData;

//------------------------------------------------------
struct sPoint {
    sPoint()
        : id(-1)
        , x(0)
        , y(0)
    {
    }
    int id;
    double x; // position of point center
    double y; // position of point center
};

//------------------------------------------------------
struct sFace {

    sFace()
        : id(-1)
        , bTypeScalar(0)
        , bTypeVelocity(0)
        , bValueX(42)
        , bValueY(42)
    {
    }
    int id;

    // grid settings
    double x;             // position of face center
    double y;             // position of face center
    double dx;            // delta x
    double dy;            // delta y
    sPoint* points[2];    // face points
    sCell* neighCells[2]; // two neighbor cells

    // boundary settings
    int bTypeScalar;      // boundary type scalar
    int bTypeVelocity;      // boundary type velocity
    double bValueX; // boundary value
    double bValueY; // boundary value

    // numerical settings
    double numFlux[2]; // numerical flux in x,y

    // physical settings
    double u;       // velocity
    double v;       // velocity
    double uNext; // velocity
    double vNext; // velocity
    
    // coefficient for simple
    double apTilde;
};

//------------------------------------------------------
struct sCell {

    sCell()
        : id(-1)
        , bTypeScalar(0)
        , bTypeVelocity(0)
        , s(0.)
    {
    }
    int id;

    // grid quantities
    double x;             // position of cell center
    double y;             // position of cell center
    sFace* faces[4];      // cell faces
    sPoint* points[4];    // cell points
    sCell* neighCells[4]; // neighbour cells

    // numeric quantities
    int bTypeScalar;      // boundary type for scalar
    int bTypeVelocity;    // boundary type for velocity
    double bValueScalar;  // boundary value
    double bValueScalarX; // boundary value
    double bValueScalarY; // boundary value
    double bValueU;       // boundary value
    double bValueV;       // boundary value

    // physical quantities
    double volume;      // cell volume
    double fluxBalance; // flux balance of cell
    double phi;         // phi at cell center
    double s;           // source term
    double p;           // pressure
    double pCorrect;    // pressure corrector
};

//------------------------------------------------------
struct sData {

    // points
    int nPoints;    // total number of points
    int nPointsX;   // number of points in x direction
    int nPointsY;   // number of points in y direction
    sPoint* points; // pointer to point array

    // cells
    int nCells;   // total number of cells
    int nCellsX;  // number of cells in x direction
    int nCellsY;  // number of cells in y direction
    sCell* cells; // pointer to cell array

    // faces
    int nFaces;   // total number of faces
    sFace* faces; // pointer to face array

    // numerical settings
    double maxTime;   // maximum physical time
    int solverMethod; // solver method
    int solverType;   // solver type
    int maxIter;      // maximum iterations
    double residuum;  // maximum residuum

    // physical settings
    double alpha;   // diffusion coefficient
    double eta;     // dynamic viscosity
    double rho;     // density
    double u;       // velocity
    double v;       // velocity
    double initPhi; // initial phi value
    double xMin;
    double xMax;
    double yMin;
    double yMax;
};

#endif
