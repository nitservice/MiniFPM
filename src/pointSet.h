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
 *   Main method for particle generation, adaptation.
 *   The initial geometry is read from a gmsh file. The boundaries are provided
 *   by "physical boundaries" in gmsh. The following physical ids are possible:
 *    1-99    for fixed boundaries (always present in simulation)
 *   100-199  for initially inactive boundaries (not filled initially - will be
 *            active by free surface contact)
 *   200-299  free surface boundaries
 *
 *************************************************************************/


#ifndef POINTSET_H
#define POINTSET_H


#include <iostream>
#include <vector>

#include "Array.h"
#include "points.h"
#include "mfpmException.h"
#include "mfpmParameters.h"
#include "refinement.h"
#include "vectorField.h"

using namespace std;

class scalarField;
class vectorField;


class pointSet
{
protected:
	struct mainGeometry
	{
		point* P;
		point* dataPoints;
		int**   line;
		double* normalX;
		double* normalY;
		bool*   active;
		int     nbLines;
		int     nbPoints;
		int     nbDataPoints;
		bool    init;
	};

	const static int activeFixedWall = 99;		// 1-99 fixed wall (active at beginning)
	const static int inactiveFixedWall = 199;   // 100-199 fixed wall (inactive at beginning)
	const static int freeSurface = 299;			// 200-299 free surface
	static mfpmParameters PS_parameter;
	static int meshCounter;
	int meshNumber;

	// global data
	vector<point*> pData;
	int NP;
	int NbBndPoints;
	double SL;
	double SLmin;
	refinement* SLrefiner;

	vector<scalarField*> sF;

	// define sector for neighbour search
	double Bbox[4];
	int nbXquads;
	int nbYquads;
	vector<int>* quadID;

	// data to store "real" neighbours
	vector<int>* neighbours;

	// data for boundary lines (defines the polygon)
	mainGeometry geometry;
	vector<int>* bndQuadID;

private:
	void updateBoundingBox();
	void initFixedBoundaryPoints(vector<point>& dataPoints);
	void parameterRepresent(int lineSeg, const point& P, double& lambda, double& dist) const;
	void findCorrespondingSegment(const point& P, int& seg, point& intersect, double eps=0.5) const;

public:
	void addGeometricPoint(int pn);
	void addFreeSurfacePoint(const point& P);
	void addInnerPoint(const point& pnt);

	void removeInnerPoint(int i);
	void removeBndPoint(int i);
	void removePoints(bool removeIt[]);
	void removeGeometricPoint(int pn);

	void movePoints(vectorField& V);
	void movePoints(const Array& Vx, const Array& Vy);

	void copyGeometry(const mainGeometry& geo);


	string name;
	void fillQuads();
	void fillBndQuads();

	pointSet(string fileName, double smoothingLength = -1.0, int nbFillSteps = -1);
	pointSet(const pointSet& ps);
	pointSet();
	~pointSet();
	void clear();
	void clearInner();

	// Functions for adaptivity
	int adaptAddInnerPoints();
	void adaptSmoothingLength();
	int adaptWallBoundaryPoints();
	int adaptFreeSurfaceBoundaryPoints();
	int adaptRemoveInnerPoints(bool outerPointsOnly = false);
	int adaptAll(bool adaptFree = false, bool adaptWall = false);
	void updateNormals();
	int updateNeighbours(int minNbP = 16, int maxNbP = 35);
	int updateNeighboursRemoveBad(int minNbP = 16, int maxNbP = 35);

	void getBoundingBox(double& xmin, double& xmax,
	                    double& ymin, double& ymax ) const;

	// "Quad" functions
	const vector<int>& getQuad(int ix, int iy, bool innerQuad = true) const;
	void getQuad(const point& P, int& ix, int& iy) const;

	// "Find neighbours" functions
	void findPossibleNeighbours(const point& P, vector<int>& Neighb, bool innerQuad=true) const;
	void findPossibleNeighbours(int i, vector<int>& Neighb, bool innerQuad=true) const;
	void findNeighbours(const point& P, vector<int>& Neighb) const;

	int  findNearestBndPoint(const point& P, bool checkInner = false) const;
	int  findNearestBndPoint(int i, bool checkInner = false) const;
	int  findNearestInnerPoint(const point& P) const;
	int  findNearestInnerPoint(int i) const;
	void findLeftRightBndPoint(int pnt, int& N1, int& N2) const;

	// Calculate coordinates w.r.t. to normal and tangent vector
	void getConeParameters(const point& P, int i, double& tau1, double& tau2) const;
	void getConeParameters(const point& P, const point& Pi, double& tau1, double& tau2) const;

	// Point in "geometric" domain (does not take free bnds into account!)
	bool pointInDomain(const point& P, double eps=0.0) const;

	// Remaining tools
	double dist(int p1, int p2) const;

	// get and set routines
	void   setBND(int i, int val);
	int    getNIDX(int i) const;
	int    getIDX(int i, int j) const;
	int    getBND(int i) const;
	Array  getBND() const;
	double getSL(const point& P) const;
	double getSL(int i) const;
	double getSL() const;
	double getMaxSL() const;
	double getMinSL() const;

	int    getMaxPoints() const;
	double getNormal(int direc, int pnt) const;
	double& getNormal(int comp, int pnt);
	Array  getNormal(int direct) const;
	void   setNormal(int direct, int pnt, double val);

	int differentBoundaries();

	void setRefiner(refinement& SF);
	void clearRefiner();

	double& x(int i);
	double& y(int i);
	double x(int i) const;
	double y(int i) const;
	Array x() const;
	Array y() const;
	const point& P(int i) const;
	point& P(int i);
	int nbPoints() const;
	int getNbBndPoints() const;

	void operator=(const pointSet& ps);
	void copyMeshWithUpdate(const pointSet& ps);
	void copyMesh(const pointSet& ps);
	void copyPointSet(const pointSet& ps, bool cG = true);

	// Add / remove pointer to scalarFields refering to this point set
	void addScalarField(scalarField* sf);
	void removeScalarField(scalarField* sf);

	void printParameters();
	void printGeoInfo();
	void printPointData();

	void write_ADOLC_Array(
	    int bnd_block_size,
	    int inner_block_size,
	    int nb_bnd_blocks,
	    int max_nb_neighbours,
	    int* bnd_block_list,
	    double* X_arr,
	    double* Y_arr,
	    int*    Plist
	);

};







#endif
