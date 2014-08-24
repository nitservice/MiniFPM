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
 *   Defines a scalarField which is an extented array depending on a pointSet
 *
 *************************************************************************/


#ifndef SCALARFIELD_H
#define SCALARFIELD_H


#include <iostream>
#include "Array.h"

using namespace std;

class pointSet;
class point;

class scalarField
{
protected:
	pointSet* mesh;
	Array data;

public:
	string name;

	scalarField();
	scalarField(pointSet &mesh);
	scalarField(pointSet &mesh, const string& name);
	scalarField(const scalarField& sf);

	~scalarField();


	// Inline
	pointSet* getMesh() const  {
		return mesh;
	};
	Array& getData() {
		return data;
	};
	Array& operator()() {
		return data;
	};
	Array operator()() const {
		return data;
	};


	double value(int i) const {
		return data.value(i);
	};
	int    getSize() const {
		return data.getSize();
	};
	double interpolate(const point& P, bool inner = true) const;
	double interpolateLS(const point& P, bool forceLowOrder = false) const;
	void changeMesh(const pointSet& meshNew);
	void assignBndValue(const scalarField& sf, int bnd);


	// Operators
	double& operator()(int i);
	double operator()(int i) const;
	double operator()(double x, double y) const;
	void operator=(const scalarField& sf);
	void operator=(const Array& v);
	void operator=(const double v);


	// Remaining stuff
	void print() const;

	void removePoint(int i);
	void addPoint(int i, double val);
	void addPoint(int i);
	void addBndPoint(int pos, const point& P);
	void addInnerPoint(const point& P);

	void printData();
};


#endif
