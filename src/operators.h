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
 *   Calculates the discretization matrix of the mathematical operators.
 *
 *************************************************************************/


#ifndef OPERATORS_H
#define OPERATORS_H

#include "mfpmMatrix.h"

class pointSet;
class scalarField;
class vectorField;
class Array;



mfpmMatrix laplace(const pointSet& mesh, bool derivBnd=false);

mfpmMatrix laplace4(const pointSet& mesh, bool derivBnd=false);

mfpmMatrix biLaplace(const pointSet& mesh);

mfpmMatrix laplace(scalarField& sf);

mfpmMatrix divEtaGrad(const scalarField &mu);

mfpmMatrix ddt(const pointSet& mesh, double tau, scalarField& oldTime);

void gradient(const pointSet& mesh, mfpmMatrix& Gx, mfpmMatrix& Gy, int order = 2);

mfpmMatrix gradNormal(const pointSet& mesh);

mfpmMatrix isource(const pointSet& mesh);

void explDerivatives(const scalarField& f, vectorField* grad, scalarField* lapF, scalarField* biLapF, bool inner=true);

scalarField smoothField(scalarField &sf, double alpha = 1.0, double smoothingLength = -1.0);

scalarField smoothField_full(scalarField &sf, double alpha = 1.0, double smoothingLength = -1.0);

#endif
