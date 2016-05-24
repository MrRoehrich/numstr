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
#include <iostream>
#include <fstream>
#include <string>

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "input.h"
#include "data.h"

//------------------------------------------------------
bool input(const char* cfgFilePath, sData* data)
{
  char line[256]=" ", token[16]=" ";
  std::cout << "\nInput:\n-------\n";

  //////////////////////
  // READ CFG FILE    //
  //////////////////////

  // open input file
  std::ifstream cfgFile(cfgFilePath);
  if(!cfgFile) { return false; }

  // read input file line by line
  while(!cfgFile.eof()) {
    cfgFile.getline(line,255);
    if(!sscanf(line,"%15s",token)){ continue; };

    if(!strcmp(token,"#")) {
      // skip comment lines
        } else if(!strcmp(token,"nCellsX")) {
      sscanf(line,"%15s %d",token,&data->nCellsX);
        } else if(!strcmp(token,"nCellsY")) {
      sscanf(line,"%15s %d",token,&data->nCellsY);
        } else if(!strcmp(token,"xMin")) {
      sscanf(line,"%15s %lf",token,&data->xMin);
        } else if(!strcmp(token,"xMax")) {
      sscanf(line,"%15s %lf",token,&data->xMax);
        } else if(!strcmp(token,"yMin")) {
      sscanf(line,"%15s %lf",token,&data->yMin);
        } else if(!strcmp(token,"yMax")) {
      sscanf(line,"%15s %lf",token,&data->yMax);
      // Add new keywords here! 
        } else {
      std::cout << "unknown token: " << token << std::endl;
    }
  }
  cfgFile.close();

  return true;
}


