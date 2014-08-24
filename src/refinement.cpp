#include "refinement.h"
#include <cmath>
#include <iostream>
#include "mfpmException.h"
#include "points.h"
#include "writeDebugVTK.h"

#include "pointSet.h"
#include "scalarField.h"


refinement::refinement(const pointSet& mesh, double factor)
{
	h = 0.5*factor*mesh.getMinSL();
	refValue = NULL;
	setNewMesh(mesh);
	double SLmax = mesh.getMaxSL();
	for (int i=0; i<nbXquads*nbYquads; i++)
		refValue[i] = SLmax;
}

refinement::~refinement()
{
	delete[] refValue;
	refValue = NULL;
}


void refinement::setNewMesh(const pointSet& mesh)
{
	mesh.getBoundingBox(xmin, xmax, ymin, ymax);
	nbXquads = int( (xmax-xmin)/h ) +2;
	nbYquads = int( (ymax-ymin)/h ) +2;
	delete[] refValue;
	refValue = new double[nbXquads*nbYquads];
}



void refinement::adapt(const scalarField& sf)
{
	setNewMesh(*sf.getMesh());
	point P;
	double SLmax = sf.getMesh()->getMaxSL();
	double SLmin = sf.getMesh()->getMinSL();

	for (int i=0; i<nbXquads; i++)
	{
		P.X = i*h + xmin;
		for (int j=0; j<nbYquads; j++)
		{
			P.Y = j*h + ymin;
			double val = max(0.0,min(abs(sf(P.X, P.Y)),1.0));
			operator()(i,j) = (1-val)*SLmax + val*SLmin;
		}
	}

}




double& refinement::operator()(int xi, int yi)
{
	int idx = xi + nbXquads*yi;
	if (idx>nbYquads*nbXquads || idx <0)
		throw mfpmExcept(2);
	return refValue[idx];
}


double& refinement::operator()(double x, double y)
{
	int ix = max(min(int( (x-xmin) / h ),nbXquads-1),0);
	int iy = max(min(int( (y-ymin) / h ),nbYquads-1),0);

	return operator()(ix, iy);
}



void refinement::writeVTK()
{
	double* X = new double[nbXquads*nbYquads];
	double* Y = new double[nbXquads*nbYquads];
	double* val = new double[nbXquads*nbYquads];

	for (int i=0; i<nbXquads; i++)
		for (int j=0; j<nbYquads; j++)
		{
			X[i+j*nbXquads] = i*h + xmin;
			Y[i+j*nbXquads] = j*h + ymin;
			val[i+j*nbXquads] = operator()(i,j);
		}
	writeDebugVTK(nbXquads*nbYquads, X, Y, val);
	delete[] X;
	delete[] Y;
	delete[] val;

}

