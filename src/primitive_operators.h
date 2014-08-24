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
 *   Primitive operators. Calculates least-squares approximation of
 *   mathematical operators at a given point.
 *
 *************************************************************************/


#ifndef PRIMITIVE_OPERATORS
#define PRIMITIVE_OPERATORS

#include <vector>

class pointSet;
class Array;
class scalarField;

extern "C" void dgels_(char* c, int *n, int* m, int* one, double *A, int* n1, double* b, int* n2, double* work, int* workl, int* INFO);
extern "C" void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda,
                        double *b, int *ldb, int *JPVT, double *RCOND, int *RANK,
                        double *work, int *lwork, int *info);


int laplaceMatEntry(const pointSet& mesh, int pnb, Array& entries);
int partialLaplaceMatEntry(const pointSet& mesh, int pnb, Array& active, Array& entries);

int laplaceOrder4MatEntry(const pointSet& mesh, int pnb, Array& entries);
int biLaplaceMatEntry(const pointSet& mesh, int pnb, Array& entries);
int gradMatEntry(const pointSet& mesh, int pnb, Array& gradX, Array& gradY, int order=2);
int partialGradMatEntry(const pointSet& mesh, int pnb, const Array& active, Array& gradX, Array& gradY, int order=2);

int neumannMatEntry(const pointSet& mesh, int pnb, Array& entries, bool forceLowOrder=false);
int partialNeumannMatEntry(const pointSet& mesh, int pnb, const Array& active, Array& entries, bool forceLowOrder=true);


void approximatedNormals(const pointSet& mesh, std::vector<int> neighb, int pnb, double& nx, double& ny);

int interpolationConstants(const scalarField& f, int pos, double* con, bool useInner=true);

int partialApproximationMatEntry(const pointSet& mesh, int pnb, Array& active, Array& entries);


int divLambdaGradMatEntry(const pointSet& mesh, int pnb, double lambda, double DxLambda, double DyLambda, const Array& larr, Array& entries);

int approximationXYMatEntry(const pointSet& mesh, int pnb, double xp, double yp, Array& entries, int order=1);


#endif
