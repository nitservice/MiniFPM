/*************************************************************************
 *
 *  This file is part of MiniFPM
 *  Copyright (C) 2014 Jan Marburger
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Description:
 *   Writes vector- and scalarFields to a VTK file.
 *
 *************************************************************************/


#ifndef POSTRESULT_H
#define POSTRESULT_H

#include <iostream>


using namespace std;

class pointSet;
class scalarField;
class vectorField;

class postResult
{
private:
	static const int nbDataSets = 20;
	pointSet* mesh;
	scalarField* sData[nbDataSets];
	vectorField* vData[nbDataSets];
	int nbScalarData;
	int nbVectorData;

	string fileName;

	void convertDoubleChar(unsigned char out[8], double in);
	void convertIntegerChar(unsigned char out[4], int in);

public:
	postResult(pointSet& mesh);
	postResult(pointSet& mesh, string fn);
	postResult(pointSet& mesh, int suffix);

	void add(scalarField& val);
	void add(vectorField& val);

	void writeVTK(const string& pfileName, bool printMessage = true);
	void writeVTK(const char* cfileName, bool printMessage = true);
	void writeVTK(int suffix, bool printMessage = true);
	void writeVTK(bool printMessage = true);

	void writeGeometryVTK(const string& fileName, bool printMessage = true);
	void writeGeometryVTK(const char* cfileName, bool printMessage = true);
	void writeGeometryVTK(int suffix, bool printMessage = true);
	void writeGeometryVTK(bool printMessage = true);

	void writeBoundaryVTK(const string &fileName, bool printMessage = true);
	void writeBoundaryVTK(const char* cfileName, bool printMessage = true);
	void writeBoundaryVTK(int suffix, bool printMessage = true);
	void writeBoundaryVTK(bool printMessage = true);

	void setSuffix(int vtk, int bnd_vtk, int geo_vtk);
	void setSuffix(int all_vtk);


};




#endif
