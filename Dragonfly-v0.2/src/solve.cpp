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

#include <math.h>
#include <stdio.h>

#include "data.h"
#include "output.h"
#include "setup.h"
#include "solve.h"

//------------------------------------------------------
double A(double Pe, sData* data)
{
    if(data->solverType == CENTRAL)
        return 1. - ABS(Pe) / 2.;
    else if(data->solverType == UPWIND)
        return 1.;
    else if(data->solverType == HYBRID)
        return MAX(0., 1. - ABS(Pe) / 2.);
    else if(data->solverType == POWER)
        return MAX(0., powf(1. - ABS(Pe) / 10., 5.));
    else if(data->solverType == EXPONENTIAL) {
        if(ABS(Pe) < 1e-10)
            return 1.;
        return ABS(Pe) / (exp(ABS(Pe)) - 1.);
    } else {
        std::cout << "\ninvalid solver type\n" << std::endl;
        return -1;
    }
}

//------------------------------------------------------
double calcPe(sData* data, double velocity, double length)
{
    if(data->eta > EPS) {
        return data->rho * velocity / data->eta * fabs(length);
    } else {
        if(velocity > 0.)
            return MAXDOUBLE;
        else
            return -MAXDOUBLE;
    }
}

//------------------------------------------------------
bool solveCalcFlux(sData* data)
{
    std::cout << "\nCalculation:\n------------\n";
    static sCell* curCell = 0;
    static sFace* ym = 0;
    static sFace* yp = 0;
    static sFace* xm = 0;
    static sFace* xp = 0;
    double dx, dy;

    double curTime = 0.;
    int curIter = 0.;
    double deltaT = data->maxTime / data->maxIter;

    // write initial conditions
    std::cout << "\nOutput... " << curIter << "\n";
    if(!output(data, curIter)) {
        std::cout << "ERROR while data output...exiting";
        getchar();
        return 1;
    }

    while(curTime < data->maxTime && curIter < data->maxIter) {

        ++curIter;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        // calcVelocityField(data, deltaT);

        if(data->solverType == CENTRAL)
            calcFluxCentral(data);
        else if(data->solverType == UPWIND)
            calcFluxUpwind(data);

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if(curCell->bTypeScalar == 1)
                continue;
            ym = curCell->faces[YM];
            yp = curCell->faces[YP];
            xm = curCell->faces[XM];
            xp = curCell->faces[XP];

            // fluxes f* (horizontally) and g* (vertically)
            dx = ym->points[P]->x - ym->points[M]->x;
            dy = ym->points[P]->y - ym->points[M]->y;
            curCell->fluxBalance = ym->numFlux[M] * dy - ym->numFlux[P] * dx;
            dx = yp->points[P]->x - yp->points[M]->x;
            dy = yp->points[P]->y - yp->points[M]->y;
            curCell->fluxBalance += -yp->numFlux[M] * dy + yp->numFlux[P] * dx;
            dx = xm->points[P]->x - xm->points[M]->x;
            dy = xm->points[P]->y - xm->points[M]->y;
            curCell->fluxBalance += -xm->numFlux[M] * dy - xm->numFlux[P] * dx;
            dx = xp->points[P]->x - xp->points[M]->x;
            dy = xp->points[P]->y - xp->points[M]->y;
            curCell->fluxBalance += xp->numFlux[M] * dy + xp->numFlux[P] * dx;

            curCell->phi += deltaT / data->rho * curCell->s;
            curCell->phi -= deltaT / (curCell->volume * data->rho) * curCell->fluxBalance;
        }

        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            if((cId + data->nCellsX) % data->nCellsX == 0) { // linker Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[XM]->u = curCell->faces[XP]->u;
                    curCell->faces[XM]->v = curCell->faces[XP]->v;
                }
            } else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) { // rechter Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[XP]->u = curCell->faces[XM]->u;
                    curCell->faces[XP]->v = curCell->faces[XM]->v;
                }
            } else if(cId < data->nCellsX) { // unterer Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[YM]->u = curCell->faces[YP]->u;
                    curCell->faces[YM]->v = curCell->faces[YP]->v;
                }
            } else if(data->nCells - (cId + 1) < data->nCellsX) { // oberer Rand
                if(curCell->bTypeVelocity == 2) {
                    curCell->faces[YP]->u = curCell->faces[YM]->u;
                    curCell->faces[YP]->v = curCell->faces[YM]->v;
                }
            }
        }

        curTime += deltaT;
        // write output
        std::cout << "Output... " << curIter << "\n";
        if(!output(data, curIter)) {
            std::cout << "ERROR while data output...exiting";
            getchar();
            return 1;
        }
    }
    return true;
}

bool solvePe(sData* data)
{
    std::cout << "\nCalculation:\n------------\n";

    double curTime = 0.;
    int curIter = 0.;
    double deltaT = data->maxTime / data->maxIter;

    // setRiverFlowBoundariers1(data);
    // setRiverFlowBoundariers2(data);
    // setRiverFlowBoundariers3(data);
    setCouetteBoundaries(data);
    // setPoiseuilleBoundaries(data);

    // write initial conditions
    std::cout << "Output... " << curIter << "\n";
    if(!output(data, curIter)) {
        std::cout << "ERROR while data output...exiting";
        getchar();
        return 1;
    }

    while(curTime < data->maxTime && curIter < data->maxIter) {
        ++curIter;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        calcVelocityField(data, deltaT);
        // calcScalarField(data, deltaT);

        curTime += deltaT;
        // write output
        std::cout << "Output... " << curIter << "\n";
        if(!output(data, curIter)) {
            std::cout << "ERROR while data output...exiting";
            getchar();
            return 1;
        }
    }

    std::cout << "\n";
    return true;
}

bool solveSimple(sData* data)
{
    double maxRes;
    double eps = data->residuum;
    sCell* curCell;
    sFace* curFace;
    sFace* faceN;
    sFace* faceS;
    sFace* faceE;
    sFace* faceW;

    double up, un, ue, us, uw;
    double vp, vn, ve, vs, vw;
    double an, ae, aw, as, ap, apTilde, b, deltaPN, deltaPT;
    double rho, dx, dy;

    double curTime = 0.;
    int curIter = 0.;
    double deltaT = data->maxTime / data->maxIter;

    // setRigidBodyBoundaries(data);
    // setDrivenCavityBoundaries(data);
    // setGartenschlauchBoundaries(data);
    // setCouetteBoundaries(data);
    setStolperdrahtBoundaries(data);

    std::cout << "Output... " << 0 << "\n";
    if(!output(data, 0)) {
        std::cout << "ERROR while data output...exiting";
        getchar();
        return 1;
    }

    // iterate over all time steps
    while(curTime < data->maxTime && curIter < data->maxIter) {
        ++curIter;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        maxRes = MAXDOUBLE;
        // reset pressure correction to 0 for next time step
        for(int cId = 0; cId < data->nCells; cId++) {
            curCell = &data->cells[cId];
            curCell->pCorrect = 0.;
        }

        // SIMPLE: while-loop
        while(maxRes > eps) {
            maxRes = 0.;

            // compute new velocties
            for(int fId = 0; fId < data->nFaces; fId++) {

                curFace = &data->faces[fId];
                if(curFace->bTypeVelocity == DIRICHLET)
                    continue;
                else if(curFace->bTypeVelocity == NEUMANN) { // mirroring velocities
                    if(curFace->dy == 0) {
                        if(curFace->neighCells[P]->neighCells[XP] != NULL) {
                            faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                            curFace->u = faceE->u;
                            curFace->v = faceE->v;
                        } else if(curFace->neighCells[P]->neighCells[XM] != NULL) {
                            faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                            curFace->u = faceW->u;
                            curFace->v = faceW->v;
                        }
                    } else if(curFace->dx == 0) {
                        if(curFace->neighCells[P] != NULL) {
                            faceE = curFace->neighCells[P]->faces[XP];
                            curFace->u = faceE->u;
                            curFace->v = faceE->v;
                        } else if(curFace->neighCells[M] != NULL) {
                            faceW = curFace->neighCells[M]->faces[XM];
                            curFace->u = faceW->u;
                            curFace->v = faceW->v;
                        }
                    }
                } else if(curFace->bTypeVelocity == INNERCELL) {
                    deltaPN = curFace->neighCells[M]->p - curFace->neighCells[P]->p;

                    if(curFace->dy == 0) {
                        faceN = curFace->neighCells[P]->faces[YP];
                        faceS = curFace->neighCells[M]->faces[YM];
                        faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                        faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                        dx = curFace->dx;
                        dy = curFace->neighCells[M]->faces[XP]->dy;
                        deltaPT =
                            (curFace->neighCells[P]->neighCells[XM]->p + curFace->neighCells[M]->neighCells[XM]->p -
                                curFace->neighCells[M]->neighCells[XP]->p - curFace->neighCells[P]->neighCells[XP]->p) /
                            4.;

                        // v
                        vp = curFace->v;
                        vn = faceN->v;
                        ve = faceE->v;
                        vs = faceS->v;
                        vw = faceW->v;
                        ue = (curFace->neighCells[M]->faces[XP]->u + curFace->neighCells[P]->faces[XP]->u) / 2.;
                        uw = (curFace->neighCells[M]->faces[XM]->u + curFace->neighCells[P]->faces[XM]->u) / 2.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, (vn + vp) / 2., ue, (vp + vs) / 2., uw, an, ae, as,
                            aw, ap, apTilde);

                        curFace->vNext =
                            (an * vn + ae * ve + as * vs + aw * vw + ap * vp + deltaPN * curFace->dx) / apTilde;
                        curFace->apTilde = apTilde;

                        //  u
                        up = curFace->u;
                        un = faceN->u;
                        ue = faceE->u;
                        us = faceS->u;
                        uw = faceW->u;
                        vn = (curFace->v + curFace->neighCells[P]->faces[YP]->v) / 2.;
                        vs = (curFace->neighCells[M]->faces[YM]->v + curFace->v) / 2.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, vn, (ue + up) / 2., vs, (up + uw) / 2., an, ae, as,
                            aw, ap, apTilde);

                        curFace->uNext = (an * un + ae * ue + as * us + aw * uw + ap * up + deltaPT * dy) / apTilde;
                        // curFace->apTilde = apTilde;
                    } else if(curFace->dx == 0) {
                        faceE = curFace->neighCells[P]->faces[XP];
                        faceW = curFace->neighCells[M]->faces[XM];
                        faceN = curFace->neighCells[M]->neighCells[YP]->faces[XP];
                        faceS = curFace->neighCells[M]->neighCells[YM]->faces[XP];
                        dx = curFace->neighCells[M]->faces[YM]->dx;
                        dy = curFace->dy;
                        deltaPT =
                            (curFace->neighCells[P]->neighCells[YM]->p + curFace->neighCells[M]->neighCells[YM]->p -
                                curFace->neighCells[M]->neighCells[YP]->p - curFace->neighCells[P]->neighCells[YP]->p) /
                            4.;

                        // v
                        vp = curFace->v;
                        vn = faceN->v;
                        ve = faceE->v;
                        vs = faceS->v;
                        vw = faceW->v;
                        ue = (curFace->u + curFace->neighCells[P]->faces[XP]->u) / 2.;
                        uw = (curFace->neighCells[M]->faces[XM]->u + curFace->u) / 2.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, (vn + vp) / 2., ue, (vp + vs) / 2., uw, an, ae, as,
                            aw, ap, apTilde);

                        curFace->vNext = (an * vn + ae * ve + as * vs + aw * vw + ap * vp + deltaPT * dx) / apTilde;
                        // curFace->apTilde = apTilde;

                        //  u
                        up = curFace->u;
                        un = faceN->u;
                        ue = faceE->u;
                        us = faceS->u;
                        uw = faceW->u;
                        vn = (curFace->neighCells[M]->faces[YP]->v + curFace->neighCells[P]->faces[YP]->v) / 2.;
                        vs = (curFace->neighCells[M]->faces[YM]->v + curFace->neighCells[P]->faces[YM]->v) / 2.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, vn, (ue + up) / 2., vs, (up + uw) / 2., an, ae, as,
                            aw, ap, apTilde);

                        curFace->uNext =
                            (an * un + ae * ue + as * us + aw * uw + ap * up + deltaPN * curFace->dy) / apTilde;
                        curFace->apTilde = apTilde;
                    }
                }
            }

            // compute p'
            for(int cId = 0; cId < data->nCells; cId++) {
                sCell* curCell = &data->cells[cId];
                if(curCell->bTypePressure == DIRICHLET)
                    continue;
                else if(curCell->bTypePressure == INNERCELL) {
                    rho = data->rho;
                    dx = curCell->faces[YP]->dx;
                    dy = curCell->faces[XP]->dy;

                    an = rho * dx * dx / curCell->faces[YP]->apTilde;
                    ae = rho * dy * dy / curCell->faces[XP]->apTilde;
                    as = rho * dx * dx / curCell->faces[YM]->apTilde;
                    aw = rho * dy * dy / curCell->faces[XM]->apTilde;

                    // Einstroemrand
                    //if(curCell->faces[XM]->bTypeVelocity == DIRICHLET)
                    //    aw = 0;

                    b = rho * ((curCell->faces[XM]->uNext - curCell->faces[XP]->uNext) * dy +
                                  (curCell->faces[YM]->vNext - curCell->faces[YP]->vNext) * dx);
                    apTilde = ae + aw + an + as;

                    curCell->pCorrect =
                        (ae * curCell->neighCells[XP]->pCorrect + as * curCell->neighCells[YM]->pCorrect +
                            aw * curCell->neighCells[XM]->pCorrect + an * curCell->neighCells[YP]->pCorrect + b) /
                        apTilde;
                    maxRes = MAX(maxRes, ABS(b));
                }
            }
            std::cout << maxRes << std::endl;

            // compute p=p*+p'
            for(int cId = 0; cId < data->nCells; cId++) {
                curCell = &data->cells[cId];
                if(curCell->bTypePressure == INNERCELL) {
                    double omega = 0.8;                      // relaxation parameter
                    curCell->p += omega * curCell->pCorrect; // now, p is the new estimate of the pressure field
                }
            }
            for(int cId = 0; cId < data->nCells; cId++) {
                curCell = &data->cells[cId];
                if(curCell->bTypePressure == NEUMANN) {
                    if((cId + data->nCellsX) % data->nCellsX == 0) // linker Rand
                        curCell->p = curCell->neighCells[XP]->p;
                    else if((cId - (data->nCellsX - 1)) % (data->nCellsX) == 0) // rechter Rand
                        curCell->p = curCell->neighCells[XM]->p;
                    else if(data->nCells - (cId + 1) < data->nCellsX) // oberer Rand
                        curCell->p = curCell->neighCells[YM]->p;
                    else if(cId < data->nCellsX) // unterer Rand
                        curCell->p = curCell->neighCells[YP]->p;

                    if(curCell->place == 1)
                        curCell->p = curCell->neighCells[XM]->p;
                    else if(curCell->place == 2)
                        curCell->p = curCell->neighCells[YP]->p;
                    else if(curCell->place == 3)
                        curCell->p = curCell->neighCells[XP]->p;
                    else if(curCell->place == 4)
                        curCell->p = curCell->neighCells[YM]->p;
                }
            }
        }

        for(int fId = 0; fId < data->nFaces; fId++) {
            curFace = &data->faces[fId];
            if(curFace->bTypeVelocity == DIRICHLET)
                continue;
            else if(curFace->bTypeVelocity == INNERCELL) {
                curFace->v = curFace->vNext;
                curFace->u = curFace->uNext;
            } else if(curFace->bTypeVelocity == NEUMANN) { // mirroring velocities
                if(curFace->dy == 0) {
                    if(curFace->neighCells[P]->neighCells[XP] != NULL) {
                        faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                        curFace->u = faceE->u;
                        curFace->v = faceE->v;
                    } else if(curFace->neighCells[P]->neighCells[XM] != NULL) {
                        faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                        curFace->u = faceW->u;
                        curFace->v = faceW->v;
                    }
                } else if(curFace->dx == 0) {
                    if(curFace->neighCells[P] != NULL) {
                        faceE = curFace->neighCells[P]->faces[XP];
                        curFace->u = faceE->u;
                        curFace->v = faceE->v;
                    } else if(curFace->neighCells[M] != NULL) {
                        faceW = curFace->neighCells[M]->faces[XM];
                        curFace->u = faceW->u;
                        curFace->v = faceW->v;
                    }
                }
            }
        }

        // write output
        std::cout << "Output... " << curIter << "\n";
        if(!output(data, curIter)) {
            std::cout << "ERROR while data output...exiting";
            getchar();
            return 1;
        }

        curTime += deltaT;
    }
    return true;
}

bool solveSimpler(sData* data)
{
    double maxRes;
    double eps = data->residuum;
    sCell* curCell;
    sFace* curFace;
    sFace* faceN;
    sFace* faceS;
    sFace* faceE;
    sFace* faceW;

    double up, un, ue, us, uw;
    double vp, vn, ve, vs, vw;
    double an, ae, aw, as, ap, apTilde, b, deltaPN, deltaPT;
    double rho, dx, dy;

    double curTime = 0.;
    int curIter = 0.;
    double deltaT = data->maxTime / data->maxIter;

    // setRigidBodyBoundaries(data);
    // setDrivenCavityBoundaries(data);
    // setGartenschlauchBoundaries(data);
    setCouetteBoundaries(data);
    // setStolperdrahtBoundaries(data);

    std::cout << "Output... " << 0 << "\n";
    if(!output(data, 0)) {
        std::cout << "ERROR while data output...exiting";
        getchar();
        return 1;
    }

    // iterate over all time steps
    while(curTime < data->maxTime && curIter < data->maxIter) {
        ++curIter;
        if(curTime + deltaT > data->maxTime)
            deltaT = data->maxTime - curTime;

        maxRes = MAXDOUBLE;
        // SIMPLER: while-loop
        while(maxRes > eps) {
            maxRes = 0.;

            // compute pseudo-velocities
            for(int fId = 0; fId < data->nFaces; fId++) {

                curFace = &data->faces[fId];
                if(curFace->bTypeVelocity == INNERCELL) {

                    if(curFace->dy == 0) {
                        faceN = curFace->neighCells[P]->faces[YP];
                        faceS = curFace->neighCells[M]->faces[YM];
                        faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                        faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                        dx = curFace->dx;
                        dy = curFace->neighCells[M]->faces[XP]->dy;

                        // vHat
                        vp = curFace->v;
                        vn = faceN->v;
                        ve = faceE->v;
                        vs = faceS->v;
                        vw = faceW->v;
                        /*ue = (curFace->neighCells[M]->faces[XP]->u + curFace->neighCells[P]->faces[XP]->u) / 2.;
                        uw = (curFace->neighCells[M]->faces[XM]->u + curFace->neighCells[P]->faces[XM]->u) / 2.;*/
                        ue = (curFace->u + curFace->neighCells[P]->faces[XP]->u +
                                 curFace->neighCells[P]->neighCells[XP]->faces[YM]->u +
                                 curFace->neighCells[M]->faces[XP]->u) /
                            4.;
                        uw = (curFace->u + curFace->neighCells[P]->faces[XM]->u +
                                 curFace->neighCells[P]->neighCells[XM]->faces[YM]->u +
                                 curFace->neighCells[M]->faces[XM]->u) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy,
                            (vp + curFace->neighCells[P]->faces[XM]->v + vn + curFace->neighCells[P]->faces[XP]->v) /
                                4.,
                            ue,
                            (vp + curFace->neighCells[M]->faces[XP]->v + vs + curFace->neighCells[M]->faces[XM]->v) /
                                4.,
                            uw, an, ae, as, aw, ap, apTilde);

                        curFace->vNext = (an * vn + ae * ve + as * vs + aw * vw + ap * vp) / apTilde;
                        curFace->apTilde = apTilde;

                        //  uHat
                        up = curFace->u;
                        un = faceN->u;
                        ue = faceE->u;
                        us = faceS->u;
                        uw = faceW->u;
                        /*vn = (curFace->v + curFace->neighCells[P]->faces[YP]->v) / 2.;
                        vs = (curFace->neighCells[M]->faces[YM]->v + curFace->v) / 2.;*/
                        vn = (curFace->v + curFace->neighCells[P]->faces[XM]->v + curFace->neighCells[P]->faces[YP]->v +
                                 curFace->neighCells[P]->faces[XP]->v) /
                            4.;
                        vs = (curFace->v + curFace->neighCells[M]->faces[XP]->v + curFace->neighCells[M]->faces[YM]->v +
                                 curFace->neighCells[M]->faces[XM]->v) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, vn,
                            (up + curFace->neighCells[P]->faces[XP]->u + ue + curFace->neighCells[M]->faces[XP]->u) /
                                4.,
                            vs,
                            (up + curFace->neighCells[M]->faces[XM]->u + uw + curFace->neighCells[P]->faces[XM]->u) /
                                4.,
                            an, ae, as, aw, ap, apTilde);

                        curFace->uNext = (an * un + ae * ue + as * us + aw * uw + ap * up) / apTilde;

                    } else if(curFace->dx == 0) {
                        faceE = curFace->neighCells[P]->faces[XP];
                        faceW = curFace->neighCells[M]->faces[XM];
                        faceN = curFace->neighCells[M]->neighCells[YP]->faces[XP];
                        faceS = curFace->neighCells[M]->neighCells[YM]->faces[XP];
                        dx = curFace->neighCells[M]->faces[YM]->dx;
                        dy = curFace->dy;

                        // vHat
                        vp = curFace->v;
                        vn = faceN->v;
                        ve = faceE->v;
                        vs = faceS->v;
                        vw = faceW->v;
                        ue = (curFace->u + curFace->neighCells[P]->faces[YP]->u + curFace->neighCells[P]->faces[XP]->u +
                                 curFace->neighCells[P]->faces[YM]->u) /
                            4.;
                        uw = (curFace->u + curFace->neighCells[M]->faces[YM]->u + curFace->neighCells[M]->faces[XM]->u +
                                 curFace->neighCells[M]->faces[YP]->u) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy,
                            (vp + curFace->neighCells[M]->faces[YP]->v + vn + curFace->neighCells[P]->faces[YP]->v) /
                                4.,
                            ue,
                            (vp + curFace->neighCells[P]->faces[YM]->v + vs + curFace->neighCells[M]->faces[YM]->v) /
                                4.,
                            uw, an, ae, as, aw, ap, apTilde);

                        curFace->vNext = (an * vn + ae * ve + as * vs + aw * vw + ap * vp) / apTilde;

                        //  uHat
                        up = curFace->u;
                        un = faceN->u;
                        ue = faceE->u;
                        us = faceS->u;
                        uw = faceW->u;
                        vn = (curFace->v + curFace->neighCells[M]->faces[YP]->v +
                                 curFace->neighCells[M]->neighCells[YP]->faces[XP]->v +
                                 curFace->neighCells[P]->faces[YP]->v) /
                            4.;
                        vs = (curFace->v + curFace->neighCells[P]->faces[YM]->v +
                                 curFace->neighCells[P]->neighCells[YM]->faces[XM]->v +
                                 curFace->neighCells[M]->faces[YM]->v) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, vn,
                            (up + curFace->neighCells[P]->faces[YP]->u + ue + curFace->neighCells[P]->faces[YM]->u) /
                                4.,
                            vs,
                            (up + curFace->neighCells[M]->faces[YM]->u + uw + curFace->neighCells[M]->faces[YP]->u) /
                                4.,
                            an, ae, as, aw, ap, apTilde);

                        curFace->uNext = (an * un + ae * ue + as * us + aw * uw + ap * up) / apTilde;
                        curFace->apTilde = apTilde;
                    }
                }
            }

            // compute p
            for(int cId = 0; cId < data->nCells; cId++) {
                sCell* curCell = &data->cells[cId];
                if(curCell->bTypePressure == DIRICHLET)
                    continue;
                else if(curCell->bTypePressure == INNERCELL) {
                    rho = data->rho;
                    dx = curCell->faces[YP]->dx;
                    dy = curCell->faces[XP]->dy;

                    an = rho * dx * dx / curCell->faces[YP]->apTilde;
                    ae = rho * dy * dy / curCell->faces[XP]->apTilde;
                    as = rho * dx * dx / curCell->faces[YM]->apTilde;
                    aw = rho * dy * dy / curCell->faces[XM]->apTilde;

                    if(curCell->faces[XM]->bTypeVelocity == DIRICHLET)
                        aw = 0;

                    b = rho * ((curCell->faces[XM]->uNext - curCell->faces[XP]->uNext) * dy +
                                  (curCell->faces[YM]->vNext - curCell->faces[YP]->vNext) * dx);
                    apTilde = ae + aw + an + as;

                    curCell->p = (ae * curCell->neighCells[XP]->p + as * curCell->neighCells[YM]->p +
                                     aw * curCell->neighCells[XM]->p + an * curCell->neighCells[YP]->p + b) /
                        apTilde;
                    maxRes = MAX(maxRes, ABS(b));
                }
            }
            std::cout << maxRes << std::endl;

            // compute new velocties
            for(int fId = 0; fId < data->nFaces; fId++) {

                curFace = &data->faces[fId];
                if(curFace->bTypeVelocity == DIRICHLET)
                    continue;
                else if(curFace->bTypeVelocity == NEUMANN) { // mirroring velocities
                    if(curFace->dy == 0) {
                        if(curFace->neighCells[P]->neighCells[XP] != NULL) {
                            faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                            curFace->u = faceE->u;
                            curFace->v = faceE->v;
                        } else if(curFace->neighCells[P]->neighCells[XM] != NULL) {
                            faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                            curFace->u = faceW->u;
                            curFace->v = faceW->v;
                        }
                    } else if(curFace->dx == 0) {
                        if(curFace->neighCells[P] != NULL) {
                            faceE = curFace->neighCells[P]->faces[XP];
                            curFace->u = faceE->u;
                            curFace->v = faceE->v;
                        } else if(curFace->neighCells[M] != NULL) {
                            faceW = curFace->neighCells[M]->faces[XM];
                            curFace->u = faceW->u;
                            curFace->v = faceW->v;
                        }
                    }
                } else if(curFace->bTypeVelocity == INNERCELL) {

                    if(curFace->dy == 0) {
                        faceN = curFace->neighCells[P]->faces[YP];
                        faceS = curFace->neighCells[M]->faces[YM];
                        faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                        faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                        dx = curFace->dx;
                        dy = curFace->neighCells[M]->faces[XP]->dy;
                        deltaPT =
                            (curFace->neighCells[P]->neighCells[XM]->p + curFace->neighCells[M]->neighCells[XM]->p -
                                curFace->neighCells[M]->neighCells[XP]->p - curFace->neighCells[P]->neighCells[XP]->p) /
                            4.;

                        // v
                        vp = curFace->v;
                        vn = faceN->v;
                        ve = faceE->v;
                        vs = faceS->v;
                        vw = faceW->v;
                        ue = (curFace->u + curFace->neighCells[P]->faces[XP]->u +
                                 curFace->neighCells[P]->neighCells[XP]->faces[YM]->u +
                                 curFace->neighCells[M]->faces[XP]->u) /
                            4.;
                        uw = (curFace->u + curFace->neighCells[P]->faces[XM]->u +
                                 curFace->neighCells[P]->neighCells[XM]->faces[YM]->u +
                                 curFace->neighCells[M]->faces[XM]->u) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy,
                            (vp + curFace->neighCells[P]->faces[XM]->v + vn + curFace->neighCells[P]->faces[XP]->v) /
                                4.,
                            ue,
                            (vp + curFace->neighCells[M]->faces[XP]->v + vs + curFace->neighCells[M]->faces[XM]->v) /
                                4.,
                            uw, an, ae, as, aw, ap, apTilde);

                        curFace->vNext =
                            (an * vn + ae * ve + as * vs + aw * vw + ap * vp + deltaPN * curFace->dx) / apTilde;
                        curFace->apTilde = apTilde;

                        //  u
                        vn = (curFace->v + curFace->neighCells[P]->faces[XM]->v + curFace->neighCells[P]->faces[YP]->v +
                                 curFace->neighCells[P]->faces[XP]->v) /
                            4.;
                        vs = (curFace->v + curFace->neighCells[M]->faces[XP]->v + curFace->neighCells[M]->faces[YM]->v +
                                 curFace->neighCells[M]->faces[XM]->v) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, vn,
                            (up + curFace->neighCells[P]->faces[XP]->u + ue + curFace->neighCells[M]->faces[XP]->u) /
                                4.,
                            vs,
                            (up + curFace->neighCells[M]->faces[XM]->u + uw + curFace->neighCells[P]->faces[XM]->u) /
                                4.,
                            an, ae, as, aw, ap, apTilde);

                        curFace->uNext = (an * un + ae * ue + as * us + aw * uw + ap * up + deltaPT * dy) / apTilde;

                    } else if(curFace->dx == 0) {
                        faceE = curFace->neighCells[P]->faces[XP];
                        faceW = curFace->neighCells[M]->faces[XM];
                        faceN = curFace->neighCells[M]->neighCells[YP]->faces[XP];
                        faceS = curFace->neighCells[M]->neighCells[YM]->faces[XP];
                        dx = curFace->neighCells[M]->faces[YM]->dx;
                        dy = curFace->dy;
                        deltaPT =
                            (curFace->neighCells[P]->neighCells[YM]->p + curFace->neighCells[M]->neighCells[YM]->p -
                                curFace->neighCells[M]->neighCells[YP]->p - curFace->neighCells[P]->neighCells[YP]->p) /
                            4.;

                        // v
                        vp = curFace->v;
                        vn = faceN->v;
                        ve = faceE->v;
                        vs = faceS->v;
                        vw = faceW->v;
                        ue = (curFace->u + curFace->neighCells[P]->faces[YP]->u + curFace->neighCells[P]->faces[XP]->u +
                                 curFace->neighCells[P]->faces[YM]->u) /
                            4.;
                        uw = (curFace->u + curFace->neighCells[M]->faces[YM]->u + curFace->neighCells[M]->faces[XM]->u +
                                 curFace->neighCells[M]->faces[YP]->u) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy,
                            (vp + curFace->neighCells[M]->faces[YP]->v + vn + curFace->neighCells[P]->faces[YP]->v) /
                                4.,
                            ue,
                            (vp + curFace->neighCells[P]->faces[YM]->v + vs + curFace->neighCells[M]->faces[YM]->v) /
                                4.,
                            uw, an, ae, as, aw, ap, apTilde);

                        curFace->vNext = (an * vn + ae * ve + as * vs + aw * vw + ap * vp + deltaPT * dx) / apTilde;

                        //  u
                        up = curFace->u;
                        un = faceN->u;
                        ue = faceE->u;
                        us = faceS->u;
                        uw = faceW->u;
                        vn = (curFace->v + curFace->neighCells[M]->faces[YP]->v +
                                 curFace->neighCells[M]->neighCells[YP]->faces[XP]->v +
                                 curFace->neighCells[P]->faces[YP]->v) /
                            4.;
                        vs = (curFace->v + curFace->neighCells[P]->faces[YM]->v +
                                 curFace->neighCells[P]->neighCells[YM]->faces[XM]->v +
                                 curFace->neighCells[M]->faces[YM]->v) /
                            4.;

                        calcCoeff(data, data->eta, deltaT, dx, dy, vn,
                            (up + curFace->neighCells[P]->faces[YP]->u + ue + curFace->neighCells[P]->faces[YM]->u) /
                                4.,
                            vs,
                            (up + curFace->neighCells[M]->faces[YM]->u + uw + curFace->neighCells[M]->faces[YP]->u) /
                                4.,
                            an, ae, as, aw, ap, apTilde);

                        curFace->uNext =
                            (an * un + ae * ue + as * us + aw * uw + ap * up + deltaPN * curFace->dy) / apTilde;
                        curFace->apTilde = apTilde;
                    }
                }
            }

            // compute p'
            for(int cId = 0; cId < data->nCells; cId++) {
                curCell = &data->cells[cId];
                if(curCell->bTypePressure == DIRICHLET)
                    continue;
                else if(curCell->bTypePressure == INNERCELL) {
                    rho = data->rho;
                    dx = curCell->faces[YP]->dx;
                    dy = curCell->faces[XP]->dy;

                    an = rho * dx * dx / curCell->faces[YP]->apTilde;
                    ae = rho * dy * dy / curCell->faces[XP]->apTilde;
                    as = rho * dx * dx / curCell->faces[YM]->apTilde;
                    aw = rho * dy * dy / curCell->faces[XM]->apTilde;

                    b = rho * ((curCell->faces[XM]->uNext - curCell->faces[XP]->uNext) * dy +
                                  (curCell->faces[YM]->vNext - curCell->faces[YP]->vNext) * dx);
                    apTilde = ae + aw + an + as;

                    curCell->pCorrect =
                        (ae * curCell->neighCells[XP]->pCorrect + as * curCell->neighCells[YM]->pCorrect +
                            aw * curCell->neighCells[XM]->pCorrect + an * curCell->neighCells[YP]->pCorrect + b) /
                        apTilde;

                    maxRes = MAX(maxRes, ABS(b));
                }
            }
            std::cout << maxRes << std::endl;

            for(int fId = 0; fId < data->nFaces; fId++) {
                curFace = &data->faces[fId];
                if(curFace->bTypeVelocity == DIRICHLET)
                    continue;
                else if(curFace->bTypeVelocity == INNERCELL) {
                    if(curFace->dx == 0) {
                        curFace->u = curFace->uNext +
                            (curFace->neighCells[M]->pCorrect - curFace->neighCells[P]->pCorrect) * curFace->dy /
                                curFace->apTilde;
                        curFace->v = curFace->vNext;
                    } else if(curFace->dy == 0) {
                        curFace->u = curFace->uNext;
                        curFace->v = curFace->vNext +
                            (curFace->neighCells[M]->pCorrect - curFace->neighCells[P]->pCorrect) * curFace->dx /
                                curFace->apTilde;
                    }
                } else if(curFace->bTypeVelocity == NEUMANN) { // mirroring velocities
                    if(curFace->dy == 0) {
                        if(curFace->neighCells[P]->neighCells[XP] != NULL) {
                            faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
                            curFace->u = faceE->u;
                            curFace->v = faceE->v;
                        } else if(curFace->neighCells[P]->neighCells[XM] != NULL) {
                            faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
                            curFace->u = faceW->u;
                            curFace->v = faceW->v;
                        }
                    } else if(curFace->dx == 0) {
                        if(curFace->neighCells[P] != NULL) {
                            faceE = curFace->neighCells[P]->faces[XP];
                            curFace->u = faceE->u;
                            curFace->v = faceE->v;
                        } else if(curFace->neighCells[M] != NULL) {
                            faceW = curFace->neighCells[M]->faces[XM];
                            curFace->u = faceW->u;
                            curFace->v = faceW->v;
                        }
                    }
                }
            }
        }

        // write output
        std::cout << "Output... " << curIter << "\n";
        if(!output(data, curIter)) {
            std::cout << "ERROR while data output...exiting";
            getchar();
            return 1;
        }

        curTime += deltaT;
    }
    return true;
}

//------------------------------------------------------
void calcFluxCentral(sData* data)
{
    static sFace* curFace = 0;

    double conv = 0;
    double diff = 0;
    double dx, dy;
    // compute numerical flux of each face
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(curFace->bTypeScalar == 0) {

            // convective part of f*
            conv = data->rho * curFace->u * (curFace->neighCells[P]->phi + curFace->neighCells[M]->phi) / 2.;
            // diffusive part of f*
            dx = curFace->neighCells[P]->x - curFace->neighCells[M]->x;
            if(dx != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dx;
            // f* = f*(conv) + f*(diff)
            curFace->numFlux[M] = conv + diff;

            // convective part of g*
            conv = data->rho * curFace->v * (curFace->neighCells[P]->phi + curFace->neighCells[M]->phi) / 2.;
            // diffusive part of g*
            dy = curFace->neighCells[P]->y - curFace->neighCells[M]->y;
            if(dy != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dy;
            // g* = g*(conv) + g*(diff)
            curFace->numFlux[P] = conv + diff;
        } else if(curFace->bTypeScalar == 2) {
            curFace->numFlux[M] = curFace->bValueX;
            curFace->numFlux[P] = curFace->bValueY;
        }
    }
}

//------------------------------------------------------
void calcFluxUpwind(sData* data)
{
    static sFace* curFace = 0;

    double conv = 0;
    double diff = 0;
    double dx, dy;
    // compute numerical flux of each face
    for(int fId = 0; fId < data->nFaces; fId++) {
        curFace = &data->faces[fId];

        if(curFace->bTypeScalar == 0) {

            // convective part of f*
            if(curFace->u >= 0)
                conv = data->rho * curFace->u * curFace->neighCells[M]->phi;
            else
                conv = data->rho * curFace->u * curFace->neighCells[P]->phi;
            // diffusive part of f*
            dx = curFace->neighCells[P]->x - curFace->neighCells[M]->x;
            if(dx != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dx;
            // f* = f*(conv) + f*(diff)
            curFace->numFlux[M] = conv + diff;

            // convective part of g*
            if(curFace->v >= 0)
                conv = data->rho * curFace->v * curFace->neighCells[M]->phi;
            else
                conv = data->rho * curFace->v * curFace->neighCells[P]->phi;
            // diffusive part of g*
            dy = curFace->neighCells[P]->y - curFace->neighCells[M]->y;
            if(dy != 0)
                diff = -data->alpha * (curFace->neighCells[P]->phi - curFace->neighCells[M]->phi) / dy;
            // g* = g*(conv) + g*(diff)
            curFace->numFlux[P] = conv + diff;
        } else if(curFace->bTypeScalar == 2) {
            curFace->numFlux[M] = curFace->bValueX;
            curFace->numFlux[P] = curFace->bValueY;
        }
    }
}

//------------------------------------------------------
void calcScalarField(sData* data, double deltaT)
{
    sCell* curCell;
    double dx, dy;
    double phiN, phiE, phiS, phiW;
    double vn, ue, vs, uw;
    double ap, an, ae, as, aw, apTilde;

    for(int cId = 0; cId < data->nCells; cId++) {
        curCell = &data->cells[cId];
        if(curCell->bTypeScalar == 1)
            continue;

        dx = curCell->faces[YP]->dx;
        dy = curCell->faces[XP]->dy;
        vn = curCell->faces[YP]->v;
        ue = curCell->faces[XP]->u;
        vs = curCell->faces[YM]->v;
        uw = curCell->faces[XM]->u;
        calcCoeff(data, data->alpha, deltaT, dx, dy, vn, ue, vs, uw, an, ae, as, aw, ap, apTilde);
        phiN = curCell->neighCells[YP]->phi;
        phiE = curCell->neighCells[XP]->phi;
        phiS = curCell->neighCells[YM]->phi;
        phiW = curCell->neighCells[XM]->phi;

        curCell->phi = (an * phiN + ae * phiE + as * phiS + aw * phiW + ap * curCell->phi) / apTilde +
            curCell->s / data->rho * deltaT;
    }
}

//------------------------------------------------------
void calcVelocityField(sData* data, double deltaT)
{
    sFace* curFace = 0;
    double deltaPN, deltaPT;
    sFace *faceN, *faceE, *faceS, *faceW;
    double dx, dy;
    double vp, vn, ve, vs, vw;
    double up, un, ue, us, uw;
    double ap, an, ae, as, aw, apTilde;

    // compute new velocties
    for(int fId = 0; fId < data->nFaces; fId++) {

        curFace = &data->faces[fId];
        if(curFace->bTypeVelocity == DIRICHLET)
            continue;
        if(curFace->bTypeVelocity == NEUMANN) { // mirroring at the right
            if(curFace->dy == 0)
                faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
            else if(curFace->dx == 0)
                faceW = curFace->neighCells[M]->faces[XM];
            curFace->u = faceW->u;
            curFace->v = faceW->v;
            continue;
        }
        deltaPN = curFace->neighCells[M]->p - curFace->neighCells[P]->p;

        if(curFace->dy == 0) {
            faceN = curFace->neighCells[P]->faces[YP];
            faceS = curFace->neighCells[M]->faces[YM];
            faceE = curFace->neighCells[P]->neighCells[XP]->faces[YM];
            faceW = curFace->neighCells[P]->neighCells[XM]->faces[YM];
            dx = curFace->dx;
            dy = curFace->neighCells[M]->faces[XP]->dy;
            deltaPT = (curFace->neighCells[P]->neighCells[XM]->p + curFace->neighCells[M]->neighCells[XM]->p -
                          curFace->neighCells[M]->neighCells[XP]->p - curFace->neighCells[P]->neighCells[XP]->p) /
                4.;

            //  u
            up = curFace->u;
            un = faceN->u;
            ue = faceE->u;
            us = faceS->u;
            uw = faceW->u;
            vn = (curFace->v + curFace->neighCells[P]->faces[YP]->v) / 2.;
            vs = (curFace->neighCells[M]->faces[YM]->v + curFace->v) / 2.;

            calcCoeff(
                data, data->eta, deltaT, dx, dy, vn, (ue + up) / 2., vs, (up + uw) / 2., an, ae, as, aw, ap, apTilde);

            curFace->u = (an * un + ae * ue + as * us + aw * uw + ap * up + deltaPT * dy) / apTilde;

            // v
            vp = curFace->v;
            vn = faceN->v;
            ve = faceE->v;
            vs = faceS->v;
            vw = faceW->v;
            ue = (curFace->neighCells[M]->faces[XP]->u + curFace->neighCells[P]->faces[XP]->u) / 2.;
            uw = (curFace->neighCells[M]->faces[XM]->u + curFace->neighCells[P]->faces[XM]->u) / 2.;

            calcCoeff(
                data, data->eta, deltaT, dx, dy, (vn + vp) / 2., ue, (vp + vs) / 2., uw, an, ae, as, aw, ap, apTilde);

            curFace->v = (an * vn + ae * ve + as * vs + aw * vw + ap * vp + deltaPN * curFace->dx) / apTilde;
        } else if(curFace->dx == 0) {
            faceE = curFace->neighCells[P]->faces[XP];
            faceW = curFace->neighCells[M]->faces[XM];
            faceN = curFace->neighCells[M]->neighCells[YP]->faces[XP];
            faceS = curFace->neighCells[M]->neighCells[YM]->faces[XP];
            dx = curFace->neighCells[M]->faces[YM]->dx;
            dy = curFace->dy;
            deltaPT = (curFace->neighCells[P]->neighCells[YM]->p + curFace->neighCells[M]->neighCells[YM]->p -
                          curFace->neighCells[M]->neighCells[YP]->p - curFace->neighCells[P]->neighCells[YP]->p) /
                4.;

            //  u
            up = curFace->u;
            un = faceN->u;
            ue = faceE->u;
            us = faceS->u;
            uw = faceW->u;
            vn = (curFace->neighCells[M]->faces[YP]->v + curFace->neighCells[P]->faces[YP]->v) / 2.;
            vs = (curFace->neighCells[M]->faces[YM]->v + curFace->neighCells[P]->faces[YM]->v) / 2.;

            calcCoeff(
                data, data->eta, deltaT, dx, dy, vn, (ue + up) / 2., vs, (up + uw) / 2., an, ae, as, aw, ap, apTilde);

            curFace->u = (an * un + ae * ue + as * us + aw * uw + ap * up + deltaPN * curFace->dy) / apTilde;

            // v
            vp = curFace->v;
            vn = faceN->v;
            ve = faceE->v;
            vs = faceS->v;
            vw = faceW->v;
            ue = (curFace->u + curFace->neighCells[P]->faces[XP]->u) / 2.;
            uw = (curFace->neighCells[M]->faces[XM]->u + curFace->u) / 2.;

            calcCoeff(
                data, data->eta, deltaT, dx, dy, (vn + vp) / 2., ue, (vp + vs) / 2., uw, an, ae, as, aw, ap, apTilde);

            curFace->v = (an * vn + ae * ve + as * vs + aw * vw + ap * vp + deltaPT * dx) / apTilde;
        }
    }
}

void calcCoeff(sData* data,
    double diffCoef,
    double deltaT,
    double dx,
    double dy,
    double vn,
    double ue,
    double vs,
    double uw,
    double& an,
    double& ae,
    double& as,
    double& aw,
    double& ap,
    double& apTilde)
{
    double rho = data->rho;

    an = diffCoef / dy * dx * A(vn * rho * dx / diffCoef, data) + MAX(-vn * rho * dx, 0.);
    as = diffCoef / dy * dx * A(vs * rho * dx / diffCoef, data) + MAX(vs * rho * dx, 0.);
    ae = diffCoef / dx * dy * A(ue * rho * dy / diffCoef, data) + MAX(-ue * rho * dx, 0.);
    aw = diffCoef / dx * dy * A(uw * rho * dy / diffCoef, data) + MAX(uw * rho * dx, 0.);

    ap = rho * dx * dy / deltaT;
    apTilde = an + as + ae + aw + ap;
}
