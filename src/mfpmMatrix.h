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
 *   Setup and handles system matrix based on a pointSet
 *   Example: mfpmMatrix A = laplace(pointSet) == -2.0;
 *            A.setDirichlet(BND_ALL, 0.0);
 *            A.solveUmfpack(result);
 *   solves the equation laplace(u) = -2 with u=0 at the boundaries.
 *
 *************************************************************************/


#ifndef MFPMMATRIX_H
#define MFPMMATRIX_H

#include "Array.h"

class pointSet;
class scalarField;


class mfpmMatrix
{
private:
	const pointSet* mesh;
	double** data;
	Array* rhs;
	int n;

public:
	mfpmMatrix(const pointSet& mesh);
	mfpmMatrix(const mfpmMatrix& mat);
	mfpmMatrix();
	~mfpmMatrix();

	void clearBnd();
	void setDirichlet(int bnd, double val);
	void setDirichlet(int bnd, const Array& val);
	void setNeumann(int bnd, double val);
	void setNeumann(int bnd, const Array& val);
	void setRobin(int bnd, double eps, double val);
	void setRobin(int bnd, double eps, const Array& val);
	void addRhs(int i, double val);
	void setRhs(int i, double val);


	double operator()(int i, int j) const;
	void operator=(const mfpmMatrix& mat);
	// Multiplying
	Array operator*(const Array& x) const;
	Array operator*(const scalarField& x) const;
	mfpmMatrix operator*(const double d) const;
	// Add
	mfpmMatrix operator+(const Array& V) const;
	mfpmMatrix operator+(const double d) const;
	mfpmMatrix operator+(const mfpmMatrix& mat) const;
	// Sub
	mfpmMatrix operator-(const mfpmMatrix& mat) const;
	mfpmMatrix operator-(const double d) const;
	mfpmMatrix operator-(const Array& V) const;

	// Inline
	double getData(int i, int j) const {
		return data[i][j];
	};
	double** getData() const {
		return data;
	};
	Array* getRhs() const {
		return rhs;
	};
	double getRhs(int i) const {
		return rhs->value(i);
	};
	const pointSet* getMesh() const {
		return mesh;
	};
	int getSize() const {
		return n;
	};

	// Friends
	friend mfpmMatrix operator*(const double d, const mfpmMatrix& mat);
	friend mfpmMatrix operator*(const Array& d, const mfpmMatrix& mat);
	friend mfpmMatrix operator+(const double d, const mfpmMatrix& mat);
	friend mfpmMatrix operator+(const Array& d, const mfpmMatrix& mat);
	friend mfpmMatrix operator-(const double d, const mfpmMatrix& mat);
	friend mfpmMatrix operator-(const Array& d, const mfpmMatrix& mat);
	friend mfpmMatrix operator==(const mfpmMatrix& mat, const Array& V);
	friend mfpmMatrix operator==(const mfpmMatrix& mat, const double d);


	void setRow(int row, Array val);
	void setDiag(int row, double val);
	void addRow(int row, Array val);
	Array getRow(int row) const;
	void regularise(double factor = 1e-4);
	void normalise();

	void solveUMFPACK(scalarField& res);



};



#endif
