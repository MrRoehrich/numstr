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
#include <fstream>
#include <string.h>
#include <math.h>
#include <iostream>

#include "input.h"
#include "data.h"

//------------------------------------------------------
bool input(const char* cfgFilePath, sData* data, int& errLine)
{
   std::cout << "\nInput:\n-------\n";

   errLine = 0;
   int lineNo;
   char line[256] = " ";
   char token[16] = " ";

   // open input file
   std::ifstream cfgFile(cfgFilePath);
   if (!cfgFile) {
      return false;
   }

   // read input file line by line
   lineNo = 0;
   while (!cfgFile.eof()) {
      lineNo++;
      cfgFile.getline(line, 255);

      // skip empty lines
      if (sscanf(line, "%15s", token) < 1) {
         continue;
      };

      // skip comment lines
      if (!strcmp(token, "#")) {
         // do nothing
      } else if (!strcmp(token, "nX")) {
         if (sscanf(line, "%15s %d", token, &data->nX) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "nY")) {
         if (sscanf(line, "%15s %d", token, &data->nY) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "xMin")) {
         if (sscanf(line, "%15s %lf", token, &data->xMin) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "xMax")) {
         if (sscanf(line, "%15s %lf", token, &data->xMax) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "solverType")) {
         if (sscanf(line, "%15s %d", token, &data->solverType) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "maxIter")) {
         if (sscanf(line, "%15s %d", token, &data->maxIter) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "maxResidual")) {
         if (sscanf(line, "%15s %lf", token, &data->maxResidual) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "facRelax")) {
         if (sscanf(line, "%15s %lf", token, &data->facRelax) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "potentialFunc")) {
         if (sscanf(line, "%15s %d", token, &data->potentialFunc) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "uInfty")) {
         if (sscanf(line, "%15s %lf", token, &data->uInfty) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "vInfty")) {
         if (sscanf(line, "%15s %lf", token, &data->vInfty) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "uIn")) {
         if (sscanf(line, "%15s %lf", token, &data->uIn) != 2) {
            errLine = lineNo;
            return false;
         };
      } else if (!strcmp(token, "magnitude")) {
         if (sscanf(line, "%15s %lf", token, &data->magnitude) != 2) {
            errLine = lineNo;
            return false;
         };
      } else {
         std::cout << "unknown token: " << token << std::endl;
         return false;
      }
   }

   cfgFile.close();

   data->x = allocGrid1Mem(data, 0);
   data->y = allocGrid1Mem(data, 0);
   data->xi = allocGrid1Mem(data, 0);
   data->et = allocGrid1Mem(data, 0);
   data->dxidx = allocGrid1Mem(data, 0);
   data->detdx = allocGrid1Mem(data, 0);
   data->dxidy = allocGrid1Mem(data, 0);
   data->detdy = allocGrid1Mem(data, 0);
   data->ddxidx = allocGrid1Mem(data, 0);
   data->ddetdx = allocGrid1Mem(data, 0);
   data->ddxidy = allocGrid1Mem(data, 0);
   data->ddetdy = allocGrid1Mem(data, 0);
   data->a1 = allocGrid1Mem(data, 0);
   data->a2 = allocGrid1Mem(data, 0);
   data->a3 = allocGrid1Mem(data, 0);
   data->a4 = allocGrid1Mem(data, 0);
   data->a5 = allocGrid1Mem(data, 0);
   data->s1 = allocGrid1Mem(data, 4711.);
   data->s2 = allocGrid1Mem(data, 4711.);
   data->ds1dxi = allocGrid1Mem(data, 0);
   data->ds1det = allocGrid1Mem(data, 0);
   data->u = allocGrid1Mem(data, 0);
   data->v = allocGrid1Mem(data, 0);
   data->cp = allocGrid1Mem(data, 0);

   return true;
}
