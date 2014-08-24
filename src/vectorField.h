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
 *   Defines a vectorField which is an extented scalarField
 *
 *************************************************************************/

#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <iostream>

using namespace std;

class scalarField;
class pointSet;

class vectorField
{
private:
	scalarField* Xval;
	scalarField* Yval;


public:
	string name;

	vectorField();
	vectorField(pointSet &mesh);
	vectorField(pointSet &mesh, const string& name);
	~vectorField();

	// Inline
	scalarField& Xcomp() {
		return *Xval;
	};
	scalarField& Ycomp() {
		return *Yval;
	};

	double& Ycomp(int i);
	double& Xcomp(int i);


};



#endif
