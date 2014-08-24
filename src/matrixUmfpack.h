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
 *   Solves A * x = rhs using UMFPACK
 *   The method 3x3 solves ( A11 A12 A13 )   (x1)    (rhs1)
 *                         ( A21 A22 A23 ) * (x2)  = (rhs2)
 *                         ( A31 A32 A33 )   (x3)    (rhs3)
 *
 *************************************************************************/


#ifndef MATRIXUMFPACK_H
#define MATRIXUMFPACK_H

#include "Array.h"

class scalarField;
class mfpmMatrix;

class matrixUmfpack
{
protected:
	int* col;
	int* row;
	double* val;
	double* rhs;

	int n;
	int nz;
public:
	matrixUmfpack(int size, int nbNonZero);
	matrixUmfpack();
	~matrixUmfpack();

	void setMatrix(const mfpmMatrix &mat);
	void setMatrix2x2(mfpmMatrix& M11, mfpmMatrix& M12,
	                  mfpmMatrix& M21, mfpmMatrix& M22 );
	void setMatrix3x3(mfpmMatrix& M11, mfpmMatrix& M12, mfpmMatrix& M13,
	                  mfpmMatrix& M21, mfpmMatrix& M22, mfpmMatrix& M23,
	                  mfpmMatrix& M31, mfpmMatrix& M32, mfpmMatrix& M33 );
	void solve(Array& Result);
	void solve(scalarField& Result);
	void solve(Array& Res1, Array& Res2);
	void solve(scalarField& Res1, scalarField& Res2);
	void solve(Array& Res1, Array& Res2, Array& Res3);
	void solve(scalarField& Res1, scalarField& Res2, scalarField& Res3);

	void operator=(const matrixUmfpack& mat);
	matrixUmfpack transpose();
	void calcDeriv
	(
	    const matrixUmfpack& matNew,
	    const matrixUmfpack& matOld,
	    double eps
	);

	double* getRhs();
	int size() {
		return n;
	};

};



#endif

