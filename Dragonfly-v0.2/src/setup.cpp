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
#include <math.h>

#include "setup.h"
#include "data.h"
#include "output.h"

//------------------------------------------------------
bool setup(sData* data)
{

   sCell* curCell = NULL;
   sFace* curFace = NULL;
   int i, j;

   /////////////////////////////////
   // construct mesh connectivity //
   /////////////////////////////////
   for (int cId = 0; cId < data->nCells; cId++) {
      // assign points to cells
      i = cId % (data->nPointsX - 1);
      j = cId / (data->nPointsX - 1);
      curCell = &data->cells[cId];
      curCell->points[XMYM] = &data->points[j * data->nPointsX + i];
      curCell->points[XPYM] = &data->points[j * data->nPointsX + i + 1];
      curCell->points[XMYP] = &data->points[(j + 1) * data->nPointsX + i];
      curCell->points[XPYP] = &data->points[(j + 1) * data->nPointsX + i + 1];

      // assign faces to cells
      curCell->faces[YM] = &data->faces[cId];
      curCell->faces[YP] = &data->faces[cId + data->nPointsX - 1];
      curCell->faces[XM] = &data->faces[(data->nPointsX - 1) * data->nPointsY + cId + j];
      curCell->faces[XP] = &data->faces[(data->nPointsX - 1) * data->nPointsY + cId + j + 1];

      curCell->faces[YM]->id = cId;
      curCell->faces[YP]->id = cId + data->nCellsX;
      curCell->faces[XM]->id = cId + data->nCellsX * (data->nCellsY + 1);
      curCell->faces[XP]->id = cId + data->nCellsX * (data->nCellsY + 1) + 1;

      // assign cells to faces
      curCell->faces[YM]->neighCells[P] = curCell;
      curCell->faces[YP]->neighCells[M] = curCell;
      curCell->faces[XM]->neighCells[P] = curCell;
      curCell->faces[XP]->neighCells[M] = curCell;

      // assign points to faces
      curCell->faces[YM]->points[M] = curCell->points[XMYM];
      curCell->faces[YM]->points[P] = curCell->points[XPYM];
      curCell->faces[YP]->points[M] = curCell->points[XMYP];
      curCell->faces[YP]->points[P] = curCell->points[XPYP];
      curCell->faces[XM]->points[M] = curCell->points[XMYM];
      curCell->faces[XM]->points[P] = curCell->points[XMYP];
      curCell->faces[XP]->points[M] = curCell->points[XPYM];
      curCell->faces[XP]->points[P] = curCell->points[XPYP];

      // assign neighboring cells to cells
      if (i != 0) {
         curCell->neighCells[XM] = &data->cells[cId - 1];
      }
      if (i != data->nPointsX - 2) {
         curCell->neighCells[XP] = &data->cells[cId + 1];
      }
      if (j != 0) {
         curCell->neighCells[YM] = &data->cells[cId - (data->nPointsX - 1)];
      }
      if (j != data->nPointsY - 2) {
         curCell->neighCells[YP] = &data->cells[cId + (data->nPointsX - 1)];
      }
   }

   /////////////////////////////////////
   // compute face centers and deltas //
   /////////////////////////////////////
   for (int fId = 0; fId < data->nFaces; fId++) {
      curFace = &data->faces[fId];
      curFace->dx = curFace->points[P]->x - curFace->points[M]->x;
      curFace->dy = curFace->points[P]->y - curFace->points[M]->y;
      curFace->x = (curFace->points[P]->x + curFace->points[M]->x) / 2.;
      curFace->y = (curFace->points[P]->y + curFace->points[M]->y) / 2.;
   }

   //////////////////////////////////////
   // compute cell centers and volumes //
   //////////////////////////////////////
   for (int cId = 0; cId < data->nCells; cId++) {
      curCell = &data->cells[cId];
      curCell->x =
          curCell->points[XMYM]->x + curCell->points[XPYM]->x + curCell->points[XMYP]->x + curCell->points[XPYP]->x;
      curCell->y /= 4.;
      curCell->y =
          curCell->points[XMYM]->y + curCell->points[XPYM]->y + curCell->points[XMYP]->y + curCell->points[XPYP]->y;
      curCell->y /= 4.;

      curCell->volume = curCell->faces[XP]->dx * curCell->faces[YP]->dy;
   }

   ///////////////////////
   // set face velocity //
   ///////////////////////
   for (int fId = 0; fId < data->nFaces; fId++) {
      curFace = &data->faces[fId];
      curFace->u = data->u;
      curFace->v = data->v;
   }

   /////////////////////////////
   // set boundary conditions //
   /////////////////////////////
   for (int cId = 0; cId < data->nCells; cId++) {
      curCell = &data->cells[cId];
      if(curCell->bType==1)
         curCell->phi = curCell->bValue;
   }


   return true;
}
