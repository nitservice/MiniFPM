#include "domainIntegration.h"

#include <cstdio>

#include "scalarField.h"
#include "pointSet.h"
#include "points.h"


#define DONT_USE_STACK


double bndIntegration(const pointSet& mesh, const scalarField& sf, int bnd)
{
	double sum = 0.0;
	int n = mesh.getNbBndPoints();
	int n1, n2;
	double dx;
	for (int i=0; i<n; i++)
	{
		if (mesh.getBND(i) == bnd)
		{
			mesh.findLeftRightBndPoint(i, n1, n2);
			// Check first neighbour
			if (n1 > i)
			{
				dx = mesh.dist(i,n1);
				sum += 0.5*(sf(i)+sf(n1)) * dx;
			}
			if (n2 > i)
			{
				dx = mesh.dist(i,n2);
				sum += 0.5*(sf(i)+sf(n2)) * dx;
			}
		}
	}
	return sum;
}



double bgGridIntegration(const pointSet& mesh, const scalarField& data)
{
	//cout << "Integrating over domain" << endl;
	// Get bounding box
	double xmin, ymin, xmax, ymax;
	mesh.getBoundingBox(xmin, xmax, ymin, ymax);

	// Get smoothing length
	double sl = mesh.getMinSL();
	double h  = 0.15*sl;

	// Calc dimensions
	int nx = (int) ((xmax-xmin)/h)+2;
	int ny = (int) ((ymax-ymin)/h)+2;


	// Build background grid
#ifndef USE_STACK
	vector<point> P;
	P.reserve(nx*ny);
#else
	point P[nx*ny];
#endif
	double val;
	for (int i=0; i<nx; i++)
	{
		for (int j=0; j<ny; j++)
		{
#ifndef USE_STACK
			point pTmp;
			pTmp.X = xmin + i*h;
			pTmp.Y = ymin + j*h;
			P.push_back(pTmp);
#else
			P[i*ny+j].X = xmin + i*h;
			P[i*ny+j].Y = ymin + j*h;
#endif
		}
	}


	// Interpolate scalarField values onto bg grid
	int n = nx*ny;
	int ix, iy;
	double sum = 0;
	for (int i=0; i<n; i++)
	{
		mesh.getQuad(P[i], ix, iy);
		const vector<int>& bndquad = mesh.getQuad(ix, iy, false);
		int bqSize = bndquad.size();
		// Point definitely in / out domain?
		if (bqSize == 0)
		{
			const vector<int>& quad = mesh.getQuad(ix, iy);
			if (quad.size() != 0)
			{
				sum += data.interpolateLS(P[i]);
			}
		}
		// No? Check in/out with normals
		else
		{
			// Find nearest bnd point
			double minDist = 1e100;
			int minPos = -1;
			for (int j=0; j<bqSize; j++)
			{
				if (minDist > mesh.P(bndquad[j]).dist(P[i]))
				{
					minDist = mesh.P(bndquad[j]).dist(P[i]);
					minPos = bndquad[j];
				}
			}
			double vecX = P[i].X - mesh.P(minPos).X;
			double vecY = P[i].Y - mesh.P(minPos).Y;

			double nDir = mesh.P(minPos).Nx*vecX + mesh.P(minPos).Ny*vecY;
			if (nDir <= 0)
			{
				sum += data.interpolateLS(P[i]);
			}
		}
	}

	sum *= h*h;

	return sum;

}
