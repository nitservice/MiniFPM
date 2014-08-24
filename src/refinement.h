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
 *   Refinement class for adaptive refinement of the point cloud.
 *
 *************************************************************************/


#ifndef REFINEMENT_H
#define REFINEMENT_H

#include "Array.h"

class pointSet;
class scalarField;

class refinement
{
protected:
	double* refValue;
	int nbXquads;
	int nbYquads;
	double h;
	double xmin, xmax;
	double ymin, ymax;

public:
	refinement(const pointSet& mesh, double factor = 1.0);

	~refinement();

	void setNewMesh(const pointSet& mesh);

	void adapt(const scalarField& sf);
	double& operator()(int xi, int yi);
	double& operator()(double x, double y);

	void writeVTK();

};



#endif
