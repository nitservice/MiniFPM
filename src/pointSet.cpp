#include "pointSet.h"
#include <fstream>
#include <cstdio>
#include <cmath>

#include "scalarField.h"
#include "primitive_operators.h"


mfpmParameters pointSet::PS_parameter;
int pointSet::meshCounter = 0;
const double slFactorForQuads = 1.0;

pointSet::pointSet(string fileName, double smoothingLength, int nbFillSteps)
{
	meshNumber = meshCounter++;
	if (smoothingLength > 0)
		SL = smoothingLength;
	else
		SL = PS_parameter.SmoothingLength;
	SLmin = PS_parameter.SmoothingLengthMin;
	nbXquads = 0;
	nbYquads = 0;
	SLrefiner  = NULL;

	ifstream gmsh(fileName.c_str());

	char ctmp[100];
	int itmp, nbElements;
	double dtmp;
	double gmshVersion;

	cout << "Reading " << fileName << "..." << endl;

	// Read $MeshFormat
	gmsh >> ctmp;
	gmsh >> gmshVersion;
	itmp = gmshVersion;
	if (itmp!=2)
	{
		throw mfpmExcept(100);
	}
	gmsh >> itmp;
	gmsh >> itmp;
	gmsh >> ctmp;
	gmsh >> ctmp;
	// Read number of nodes
	gmsh >> NP;
	cout << " -> Nb. of points: " << NP << endl;
	pData.reserve(NP);

	geometry.P = new point[NP];
	for (int i=0; i<NP; i++)
	{
		gmsh >> itmp;
		gmsh >> geometry.P[i].X;
		gmsh >> geometry.P[i].Y;
		gmsh >> dtmp;
	}
	int oldNP = NP;
	// Read nb of elements
	gmsh >> ctmp;
	gmsh >> ctmp;
	gmsh >> nbElements;

	// Read boundary elements
	geometry.active    = new bool[nbElements];
	geometry.line      = new int*[nbElements];
	geometry.normalX   = new double[nbElements];
	geometry.normalY   = new double[nbElements];
	geometry.nbLines   = nbElements;
	geometry.nbPoints  = NP;
	for (int i=0; i<nbElements; i++)
		geometry.line[i] = new int[3];
	for (int i=0; i<nbElements; i++)
	{
		gmsh >> itmp;
		gmsh >> itmp;
		if (itmp != 1)
			throw mfpmExcept(101);
		gmsh >> itmp;
		gmsh >> geometry.line[i][0];
		gmsh >> itmp;
		if (gmshVersion < 2.2)
			gmsh >> itmp;
		gmsh >> itmp;
		geometry.line[i][1] = itmp - 1;
		geometry.P[itmp-1].id = i; //geometry.line[i][0];
		gmsh >> itmp;
		geometry.line[i][2] = itmp - 1;
		if ( geometry.line[i][0] > activeFixedWall
		        && geometry.line[i][0] <= inactiveFixedWall )
			geometry.active[i] = false;
		else
			geometry.active[i] = true;
	}

	// Calculate geometry normals
	double normalLength;
	for (int i=0; i<nbElements; i++)
	{
		geometry.normalX[i] = geometry.P[geometry.line[i][1]].Y - geometry.P[geometry.line[i][2]].Y;
		geometry.normalY[i] = geometry.P[geometry.line[i][2]].X - geometry.P[geometry.line[i][1]].X;
		normalLength = sqrt( pow(geometry.normalX[i],2) + pow(geometry.normalY[i],2) );
		geometry.normalX[i] /= -normalLength;
		geometry.normalY[i] /= -normalLength;
	}

	// Fill present mfpm points
	vector<point> dataPoints;
	point ptmp(geometry.P[geometry.line[0][1]]);
	ptmp.id = 0;
	ptmp.Nx  = geometry.normalX[0];
	ptmp.Ny  = geometry.normalY[0];
	ptmp.Bnd = geometry.line[0][0];
	ptmp.meshID = meshNumber;
	dataPoints.push_back(ptmp);

	NP = 0;
	nbXquads = 0;
	nbYquads = 0;
	NbBndPoints = NP;
	quadID = NULL;
	bndQuadID = NULL;
	neighbours = NULL;
	geometry.init = true;

	initFixedBoundaryPoints(dataPoints);
	int cnt = 0;
	for (int i=0; i<dataPoints.size(); i++)
		if (dataPoints[i].Bnd <= inactiveFixedWall)
			cnt++;
	geometry.dataPoints = new point[cnt];
	geometry.nbDataPoints = cnt;
	// Finally add points
	cnt = 0;
	for (int i=0; i<dataPoints.size(); i++)
	{
		if (dataPoints[i].Bnd <= inactiveFixedWall)
		{
			geometry.dataPoints[cnt] = dataPoints[i];
			geometry.dataPoints[cnt].id = cnt;
			cnt++;
		}
	}
	// Add geometry points to pointset
	for (int i=0; i<geometry.nbDataPoints; i++)
		if (geometry.dataPoints[i].Bnd < activeFixedWall)
			addGeometricPoint(i);
		else
			geometry.dataPoints[i].active = false;
	// Add free surface points
	for (int i=0; i<dataPoints.size(); i++)
		if (dataPoints[i].Bnd>inactiveFixedWall)
			addFreeSurfacePoint(dataPoints[i]);

	// Only if required fill inner domain
	if (nbFillSteps != 0)
	{
		// Add one seed-point
		// Find first active segment
		int seg = 0;
		for (int i=0; i<geometry.nbLines; i++)
			if (geometry.active[i])
			{
				seg = i;
				break;
			}
		point Ptmp = 0.5*(geometry.P[geometry.line[seg][1]] + geometry.P[geometry.line[seg][2]]);
		Ptmp.X -= 0.2*SL*geometry.normalX[seg];
		Ptmp.Y -= 0.2*SL*geometry.normalY[seg];
		addInnerPoint(Ptmp);
		//~
		//~ seg = geometry.nbPoints-1;
		//~ Ptmp = 0.5*(geometry.P[geometry.line[seg][1]] + geometry.P[geometry.line[seg][2]]);
		//~ Ptmp.X -= 0.2*SL*geometry.normalX[seg];
		//~ Ptmp.Y -= 0.2*SL*geometry.normalY[seg];
		//~ addInnerPoint(Ptmp);

		fillQuads();

		cout << " -> Fill domain..." << endl;
		int nPointsAdded = adaptAddInnerPoints();
		cout << "    -> Points added: " << nPointsAdded << endl;
		if (nbFillSteps < 0)
			cnt = -1000;
		else
			cnt = 1;
		while ( nPointsAdded > 0 && cnt < nbFillSteps)
		{
			nPointsAdded = adaptAddInnerPoints();
			cout << "    -> Points added: " << nPointsAdded << endl;
			cnt++;
		}
		cout << " -> done!" << endl;
	}

	// Set geometry active flag
	for (int i=0; i<nbElements; i++)
		if (geometry.line[i][0]>inactiveFixedWall)
		{
			geometry.active[i]=false;		// Set to false !
		}
		else
		{
			geometry.active[i]=true;
		}
	geometry.init = false;

	if (nbFillSteps != 0)
	{
		updateBoundingBox();
		updateNeighbours();
	}

}



pointSet::pointSet(const pointSet& ps)
{
	meshNumber = meshCounter++;
	quadID = NULL;
	bndQuadID = NULL;
	neighbours = NULL;
	geometry.nbLines = 0;
	NP = 0;
	NbBndPoints = 0;
	copyPointSet(ps);
}

pointSet::pointSet()
{
	meshNumber = meshCounter++;
	quadID = NULL;
	bndQuadID = NULL;
	neighbours = NULL;
	geometry.nbLines = 0;
	NP = 0;
	NbBndPoints = 0;
}

pointSet::~pointSet()
{
	clearInner();
	clear();
}


void pointSet::clearInner()
{
	if (NP>0)
	{
		for (int i=0; i<NP; i++)
		{
			delete pData[i];
		}
		pData.clear();
		delete[] quadID;
		delete[] neighbours;
		delete[] bndQuadID;
		quadID = NULL;
		neighbours = NULL;
		bndQuadID = NULL;
		NP = 0;
		NbBndPoints = 0;
		nbXquads = 0;
		nbYquads = 0;
	}
}

void pointSet::copyPointSet(const pointSet& ps, bool cG)
{
	for (int i=0; i<NP; i++)
		if (pData[i]->meshID != meshNumber)
			cout << "ERROR"<< endl;
	clearInner();

	NP = ps.nbPoints();
	NbBndPoints = ps.getNbBndPoints();
	SL = ps.SL;
	SLmin = ps.SLmin;
	SLrefiner = ps.SLrefiner;

	// Copy geometry
	clear();
	copyGeometry(ps.geometry);

	// Copy point data vector
	pData.clear();
	pData.reserve(NP);
	for (int i=0; i<NP; i++)
	{
		point* pTmp = new point(ps.P(i));
		pTmp->meshID = meshNumber;
		pData.push_back(pTmp);
	}
	// Copy bounding box
	for (int i=0; i<4; i++)
		Bbox[i] = ps.Bbox[i];

	// Copy quads
	nbXquads = ps.nbXquads;
	nbYquads = ps.nbYquads;

	quadID = new vector<int>[nbXquads*nbYquads];
	for (int i=0; i<nbXquads*nbYquads; i++)
	{
		quadID[i] = ps.quadID[i];
	}
	bndQuadID = new vector<int>[nbXquads*nbYquads];
	for (int i=0; i<nbXquads*nbYquads; i++)
	{
		bndQuadID[i] = ps.bndQuadID[i];
	}

	// Copy neighbours
	neighbours = new vector<int>[NP];
	for (int i=0; i<NP; i++)
	{
		neighbours[i] = ps.neighbours[i];
	}

}


void pointSet::copyGeometry(const mainGeometry& geo)
{
	if (geometry.nbLines <= 0)
	{
		geometry.nbLines      = geo.nbLines;
		geometry.nbPoints     = geo.nbPoints;
		geometry.nbDataPoints = geo.nbDataPoints;
		geometry.active       = new bool[geo.nbLines];
		geometry.line         = new int*[geo.nbLines];
		for (int i=0; i<geo.nbLines; i++)
			geometry.line[i] = new int[3];
		geometry.normalX      = new double[geo.nbLines];
		geometry.normalY      = new double[geo.nbLines];
		geometry.P            = new point[geo.nbPoints];
		geometry.dataPoints   = new point[geo.nbDataPoints];
		geometry.init         = geo.init;

	}

	for (int i=0; i<geo.nbLines; i++)
	{
		geometry.active[i] = geo.active[i];
		for (int j=0; j<3; j++)
			geometry.line[i][j] = geo.line[i][j];
		geometry.normalX[i]    = geo.normalX[i];
		geometry.normalY[i]    = geo.normalY[i];
	}

	for (int i=0; i<geo.nbPoints; i++)
	{
		geometry.P[i] = geo.P[i];
	}

	for (int i=0; i<geo.nbDataPoints; i++)
	{
		geometry.dataPoints[i] = geo.dataPoints[i];

	}

}


void pointSet::clear()
{
	if (geometry.nbLines>0)
	{
		delete[] geometry.P;
		delete[] geometry.dataPoints;
		for (int i=0; i<geometry.nbLines; i++)
			delete[] geometry.line[i];
		delete[] geometry.line;
		delete[] geometry.normalX;
		delete[] geometry.normalY;
		delete[] geometry.active;
		geometry.P = NULL;
		geometry.dataPoints = NULL;
		geometry.line = NULL;
		geometry.normalX = NULL;
		geometry.normalY = NULL;
		geometry.active = NULL;
		geometry.nbLines = 0;
		geometry.nbDataPoints = 0;
		geometry.nbPoints = 0;
	}
}








int pointSet::adaptAddInnerPoints()
{
	const double nSeg = 45;
	double r    = PS_parameter.AddInner_Radius*SL;
	point pTmp;
	int count = 0;
	vector<int> posNeighb;
	fillQuads();
	fillBndQuads();
	int startPoint;
	startPoint = NbBndPoints;
	int oldNP = NP;


	for (int i=startPoint; i<NP; i++)
	{
		//~ if (!geometry.init && NP>8.0*oldNP)
		//~ return count;

		findPossibleNeighbours(i, posNeighb);
		// Set hole radius
		r = PS_parameter.AddInner_Radius*getSL(i);
		// Fill holes
		for (int seg=0; seg<nSeg; seg++)
		{
			pTmp.X = pData[i]->X + r*sin( seg/nSeg * 2*M_PI );
			pTmp.Y = pData[i]->Y + r*cos( seg/nSeg * 2*M_PI );
			if (!pointInDomain(pTmp, PS_parameter.AddInner_DeathZone))
				continue;
			bool add = true;
			for (int j=0; j<posNeighb.size(); j++)
			{
				if ( pData[posNeighb[j]]->dist(pTmp) < r-0.0001)
				{
					add = false;
					break;
				}
			}
			if (add)
			{
				pTmp.Bnd = 0;
				addInnerPoint(pTmp);
				fillQuads();
				fillBndQuads();
				count++;
				break;
			}
		}
	}
	return count;
}


int pointSet::adaptRemoveInnerPoints(bool outerPointsOnly)
{
	fillQuads();
	fillBndQuads();
	bool remove[NP];
	bool checked[NP];
	for (int i=0; i<NP; i++)
	{
		remove[i]  = false;
		checked[i] = false;
	}

	for (int i=NbBndPoints; i<NP; i++)
	{
		if (!checked[i] && !pointInDomain(*pData[i], PS_parameter.RemInner_DeathZone))
		{
			remove[i] = true;
			checked[i] = true;
		}
		else if (!outerPointsOnly)
		{
			if (checked[i])
				continue;
			vector<int> neighb;
			findPossibleNeighbours(i, neighb);
			double minDist = 1e100;
			int minPos = -1;
			for (int j=0; j<neighb.size(); j++)
			{
				if ( neighb[j] != i && !checked[neighb[j]] && pData[i]->dist(*pData[neighb[j]]) < minDist )
				{
					minDist = pData[i]->dist(*pData[neighb[j]]);
					minPos  = neighb[j];
				}
			}
			if (minPos>=0 && minDist < PS_parameter.RemInner_Distance*getSL(i))
			{
				remove[i]       = true;
				checked[i]      = true;
				checked[minPos] = true;
			}
		}
	}

	removePoints(remove);

	int count=0;
	for (int i=0; i<NP; i++)
		if (remove[i])
			count++;
	return count;
}


void pointSet::initFixedBoundaryPoints(vector<point>& dataPoints)
{

	const double addDist = PS_parameter.InitFixed_Dist*SL;
	double h = 0;
	geometry.nbDataPoints = 0;
	for (int seg=0; seg<geometry.nbLines; seg++)
	{
		if (seg>0 && geometry.line[seg-1][2] != geometry.line[seg][1] )
			h = 0;
		int A = geometry.line[seg][1];
		int B = geometry.line[seg][2];
		double length = geometry.P[A].dist(geometry.P[B]);
		if ( length + h < addDist)
		{
			h += length;
		}
		else
		{
			double hLine = addDist - h;
			point P;
			while ( hLine < length - 0.0*addDist)
			{
				P.X = geometry.P[A].X + hLine/length*(geometry.P[B].X-geometry.P[A].X);
				P.Y = geometry.P[A].Y + hLine/length*(geometry.P[B].Y-geometry.P[A].Y);
				P.Bnd = geometry.line[seg][0];
				P.Nx = geometry.normalX[seg];
				P.Ny = geometry.normalY[seg];
				P.id = seg;
				P.active = geometry.active[seg];
				P.meshID = meshNumber;
				dataPoints.push_back(P);
				hLine += addDist;
			}
			h = length-(hLine-addDist);
		}
	}

}


int pointSet::adaptWallBoundaryPoints()
{
	fillQuads();
	fillBndQuads();
	int ix, iy;
	double tx, nx;
	for (int i=0; i<geometry.nbDataPoints; i++)
	{
		vector<int> quad;
		findPossibleNeighbours(geometry.dataPoints[i], quad);
		bool found1 = false;
		bool found2 = false;
		for (int j=0; j<quad.size(); j++)
		{
			if (i==quad[j])
				continue;
			getConeParameters(geometry.dataPoints[i],*pData[quad[j]],tx,nx);
			if (   abs(tx) < PS_parameter.WallBnd_ConeLength*SL
			        && abs(nx) < PS_parameter.WallBnd_ConeLength*SL
			        && nx      < 0
			        && abs(nx) < PS_parameter.WallBnd_ConeAngle * abs(tx)
			   )
			{
				if ( tx > 0 )
				{
					found1 = true;
				}
				if ( tx < 0)
				{
					found2 = true;
				}
			}
		}
		if ( (!found1 | !found2) && geometry.dataPoints[i].active)
		{
			removeGeometricPoint(i);
			fillQuads();
			fillBndQuads();
		}
		else if ((found1&found2) && !geometry.dataPoints[i].active)
		{
			addGeometricPoint(i);
			fillQuads();
			fillBndQuads();
		}
	}
}


int pointSet::adaptFreeSurfaceBoundaryPoints()
{
	fillQuads();
	fillBndQuads();

	// Remove particles which are too close
	bool checked[NbBndPoints];
	bool remove[NP];
	for (int i=0; i<NP; i++)
		remove[i] = false;
	for (int i=0; i<NbBndPoints; i++)
		if (pData[i]->Bnd > inactiveFixedWall)
			checked[i] = false;
		else
			checked[i] = true;
	int Nr, Nl;

	for(int i=0; i<NbBndPoints; i++)
	{
		if (pData[i]->Bnd > inactiveFixedWall)
		{
			// Point out of fixed walls?
			if (!pointInDomain(*pData[i], PS_parameter.FreeSurf_DeathZone))
			{
				checked[i] = true;
				remove[i]  = true;
				continue;
			}
			// No inner point in neighbourhood?
			else
			{
				vector<int> qPoint;
				findPossibleNeighbours(i, qPoint);
				double dDist = 1e100;
				double dist;
				int cntInnerPoints = 0;
				for (int j=0; j<qPoint.size(); j++)
				{
					dist = pData[i]->dist(*pData[qPoint[j]]);
					if (dist < PS_parameter.FreeSurf_InnerDist*getSL(i) && pData[qPoint[j]]->Bnd==0)
						cntInnerPoints++;
				}
				if (cntInnerPoints < PS_parameter.FreeSurf_InnerCnt)
				{
					checked[i] = true;
					remove[i]  = true;
					continue;
				}
			}

			findLeftRightBndPoint(i, Nr, Nl);
			checked[i] = true;
			if (!checked[Nr]
			        && pData[i]->dist(*pData[Nr]) < PS_parameter.FreeSurf_RemDist*getSL(i)
			        && Nr != i)
			{
				remove[Nr] = true;
				checked[Nr]  = true;
			}
			if (!checked[Nl]
			        && pData[i]->dist(*pData[Nl]) < PS_parameter.FreeSurf_RemDist*getSL(i)
			        && Nl != i )
			{
				remove[Nl] = true;
				checked[Nl]  = true;
			}

			//Check for merging
			vector<int> qPoint;
			findPossibleNeighbours(i, qPoint, false);
			for (int j=0; j<qPoint.size(); j++)
			{
				if (qPoint[j]!=i)
				{
					if (pData[i]->dist(*pData[qPoint[j]])<PS_parameter.FreeSurf_RemDist*getSL(i))
					{
						if (pData[i]->Nx*pData[qPoint[j]]->Nx + pData[i]->Ny*pData[qPoint[j]]->Ny < 0)
						{
							remove[i] = true;
							remove[qPoint[j]] = true;
						}
					}
				}
			}
		}
	}
	removePoints(remove);
	fillQuads();
	fillBndQuads();


	// Add new free surface points
	int oldNbBndPoints = NbBndPoints;
	int seg;
	point pTmp;
	point pTmp2;
	for (int i=0; i<oldNbBndPoints; i++)
	{
		if (pData[i]->Bnd >inactiveFixedWall )
		{
			findLeftRightBndPoint(i, Nl, Nr);
			if ( pData[i]->dist(*pData[Nr]) > PS_parameter.FreeSurf_AddDist*getSL(i))
			{
				// Other boundary point -> add midpoint
				if (pData[Nr]->Bnd > inactiveFixedWall)
				{
					pTmp = 0.5*(*pData[i] + *pData[Nr]);
					pTmp.Bnd = pData[i]->Bnd;
					pTmp.Nx  = pData[i]->Nx;
					pTmp.Ny  = pData[i]->Ny;
					addFreeSurfacePoint(pTmp);
					fillBndQuads();
				}
				// Wall contact? -> add midpoint by wall normal
				else
				{
					findCorrespondingSegment(*pData[i], seg, pTmp2, 2.0);
					if (seg<0 || pData[i]->dist(pTmp2) >= pData[i]->dist(*pData[Nr]) )
					{
						pTmp = 0.5*(*pData[i] + *pData[Nr]);
						pTmp.Bnd = pData[i]->Bnd;
						pTmp.Nx  = pData[i]->Nx;
						pTmp.Ny  = pData[i]->Ny;
						addFreeSurfacePoint(pTmp);
						fillBndQuads();
						continue;
					}
					pTmp = 0.5*(*pData[i] + pTmp2);
					pTmp.Bnd = pData[i]->Bnd;
					pTmp.Nx  = pData[i]->Nx;
					pTmp.Ny  = pData[i]->Ny;
					addFreeSurfacePoint(pTmp);
					fillBndQuads();
				}
			}
			if (  pData[i]->dist(*pData[Nl]) > PS_parameter.FreeSurf_AddDist*getSL(i))
			{
				if (pData[Nl]->Bnd > inactiveFixedWall)
				{
					pTmp = 0.5*(*pData[i] + *pData[Nl]);
					pTmp.Bnd = pData[i]->Bnd;
					pTmp.Nx  = pData[i]->Nx;
					pTmp.Ny  = pData[i]->Ny;
					addFreeSurfacePoint(pTmp);
					fillBndQuads();
				}
				// Wall contact? -> add midpoint by wall normal
				else
				{
					findCorrespondingSegment(*pData[i], seg, pTmp2, 2.0);
					if (seg<0 || pData[i]->dist(pTmp2) >= pData[i]->dist(*pData[Nr]))
					{
						pTmp = 0.5*(*pData[i] + *pData[Nl]);
						pTmp.Bnd = pData[i]->Bnd;
						pTmp.Nx  = pData[i]->Nx;
						pTmp.Ny  = pData[i]->Ny;
						addFreeSurfacePoint(pTmp);
						fillBndQuads();
						continue;
					}
					pTmp = 0.5*(*pData[i] + pTmp2);
					pTmp.Bnd = pData[i]->Bnd;
					pTmp.Nx  = pData[i]->Nx;
					pTmp.Ny  = pData[i]->Ny;
					addFreeSurfacePoint(pTmp);
					fillBndQuads();
				}
			}
		}
	}

	updateNormals();

}



int pointSet::adaptAll(bool adaptFixed, bool adaptWall)
{
	if (adaptWall)
		adaptWallBoundaryPoints();
	if (adaptFixed)
		adaptFreeSurfaceBoundaryPoints();
	int slIncr = updateNeighbours();
	cout << "High change of SL: " << slIncr << endl;
	if (slIncr > -1)
	{
		bool outerPointsOnly = false;
		if (slIncr < -3)
			outerPointsOnly = true;
		cout << "Remove points: " << flush;
		int nrem = adaptRemoveInnerPoints(outerPointsOnly);
		cout <<  nrem << flush;
		while (nrem>4)
		{
			nrem = adaptRemoveInnerPoints();
			cout << " " << nrem << flush;
		}
		cout << endl;
		cout << "Add points:    " << flush;
		int nadd = adaptAddInnerPoints();
		cout << nadd << flush;
		while (nadd>4)
		{
			nadd = adaptAddInnerPoints();
			cout << " " << nadd << flush;
		}
		cout << "  --> " << nbPoints() << endl;
		slIncr = updateNeighbours();
		cout << "High change of SL: " << slIncr << endl;
	}
	return 0;
}



void pointSet::parameterRepresent(int lineSeg, const point& P, double& lambda, double& dist) const
{
	double P1x = geometry.P[geometry.line[lineSeg][1]].X;
	double P1y = geometry.P[geometry.line[lineSeg][1]].Y;
	double P2x = geometry.P[geometry.line[lineSeg][2]].X;
	double P2y = geometry.P[geometry.line[lineSeg][2]].Y;
	double deltaX = P2x-P1x;
	double deltaY = P2y-P1y;
	lambda = 1.0/(deltaX*deltaX+deltaY*deltaY)*((deltaX*P.X+deltaY*P.Y) - (deltaX*P1x + deltaY*P1y));
	double dx = P.X - (P1x + lambda*deltaX);
	double dy = P.Y - (P1y + lambda*deltaY);
	dist = sqrt( dx*dx + dy*dy );
	if ( dx*geometry.normalX[lineSeg] + dy*geometry.normalY[lineSeg] > 0 )
		dist *= -1.0;
}


void pointSet::findCorrespondingSegment(const point& P, int& seg, point& intersect, double eps) const
{
	double lambda, dist;
	double minDist = eps*SL;
	int    rSeg = -1;
	double rLambda;
	for (int lSeg=0; lSeg<geometry.nbLines; lSeg++)
	{
		if (geometry.active[lSeg])
		{
			parameterRepresent(lSeg, P, lambda, dist);
			if (lambda>=0
			        && lambda < 1
			        && abs(dist) < minDist)
			{
				rSeg = lSeg;
				rLambda = lambda;
				minDist = abs(dist);
			}
		}
	}

	seg = rSeg;
	if (rSeg<0)
	{
		return;
	}

	int p1 = geometry.line[rSeg][1];
	int p2 = geometry.line[rSeg][2];
	intersect.X = geometry.P[p1].X
	              + rLambda*(geometry.P[p2].X - geometry.P[p1].X);
	intersect.Y = geometry.P[p1].Y
	              + rLambda*(geometry.P[p2].Y - geometry.P[p1].Y);
	intersect.Bnd = geometry.line[rSeg][0];


}


bool pointSet::pointInDomain(const point& P, double eps) const
{
	double lambda, dist;
	double minDist = 1e100;
	int    seg = -1;
	double segLambda;
	int    ix, iy;
	vector<int> bndquad;
	findPossibleNeighbours(P, bndquad, false);
	// No boundary point close to P -> point in domain
	if (bndquad.size() == 0)
		return true;
	else if (!geometry.init)
	{
		bool onlyFree = true;
		for (int i=0; i<bndquad.size(); i++)
		{
			if (pData[bndquad[i]]->Bnd <= inactiveFixedWall)
				onlyFree=false;
		}
		if (onlyFree)
		{
			return true;
		}
	}
	// else check for segment
	for (int lSeg=0; lSeg<geometry.nbLines; lSeg++)
	{
		if (geometry.active[lSeg])
		{
			parameterRepresent(lSeg, P, lambda, dist);
			if (lambda>=0 && lambda <= 1
			        && abs(minDist) > abs(dist))
			{
				seg       = lSeg;
				minDist   = dist;
				segLambda = lambda;
			}
		}
	}
	// TODO: check for seg found?
	if (geometry.init && seg < 0)
		return false;

	if ( minDist < eps*SL)
	{
		return false;
	}

	if (seg<0)
		return false;

	return true;
}



void pointSet::addGeometricPoint(int pn)
{
	if (pn<0 || pn>geometry.nbDataPoints)
		throw mfpmExcept(2);
	geometry.dataPoints[pn].active = true;
	if (!geometry.init)
	{
		fillBndQuads();
	}
	for (int i=0; i<sF.size(); i++)
	{
		sF[i]->addBndPoint(0, geometry.dataPoints[pn]);
	}
	point* pTmp = new point(geometry.dataPoints[pn]);
	pTmp->meshID = meshNumber;
	pData.insert(pData.begin(), pTmp);
	NbBndPoints++;
	NP++;

}

void pointSet::addFreeSurfacePoint(const point& P)
{
	point* pTmp = new point(P);
	pTmp->meshID = meshNumber;
	pTmp->SL     = SL;
	if(!geometry.init)
		fillBndQuads();
	for (int i=0; i<sF.size(); i++)
	{
		sF[i]->addBndPoint(NbBndPoints, *pTmp);
	}
	pData.insert(pData.begin()+NbBndPoints, pTmp);
	NbBndPoints++;
	NP++;
}



void pointSet::removeGeometricPoint(int pn)
{
	if (pn<0 || pn>geometry.nbDataPoints)
		throw mfpmExcept(2);
	if (!geometry.dataPoints[pn].active)
		throw mfpmExcept(5);
	geometry.dataPoints[pn].active = false;
	int geoPn = 0;
	while ( pn != pData[geoPn]->id )
	{
		geoPn++;
		if (geoPn>=NbBndPoints)
			throw mfpmExcept(4);
	}
	delete pData[geoPn];
	pData[geoPn] = NULL;
	pData.erase(pData.begin()+geoPn);
	for (int i=0; i<sF.size(); i++)
	{
		sF[i]->removePoint(geoPn);
	}
	NbBndPoints--;
	NP--;

}

void pointSet::getBoundingBox(
    double& xmin, double& xmax,
    double& ymin, double& ymax
) const
{
	xmin = pData[0]->X;
	xmax = pData[0]->X;
	ymin = pData[0]->Y;
	ymax = pData[0]->Y;
	// mfpm points
	for (int i=1; i<NP; i++)
	{
		double x=pData[i]->X;
		double y=pData[i]->Y;
		if (xmin>x)
			xmin = x;
		if (xmax<x)
			xmax = x;
		if (ymin>y)
			ymin = y;
		if (ymax<y)
			ymax = y;
	}

}

void pointSet::updateBoundingBox()
{
	getBoundingBox(Bbox[0], Bbox[2], Bbox[1], Bbox[3]);
}


void pointSet::fillQuads()
{
	updateBoundingBox();
	int nx = int( (Bbox[2]-Bbox[0])/(slFactorForQuads*SL) );
	int ny = int( (Bbox[3]-Bbox[1])/(slFactorForQuads*SL) );
	if (nbXquads != nx  ||  nbYquads != ny )
	{
		if (quadID!=NULL)
			delete[] quadID;
		quadID = new vector<int>[nx*ny];
		nbXquads = nx;
		nbYquads = ny;
		fillBndQuads();
	}
	else
	{
		for (int i=0; i<nx*ny; i++)
		{
			quadID[i].clear();
		}
	}
	int ix,iy;
	for (int i=0; i<NP; i++)
	{
		getQuad(*pData[i], ix, iy);
		quadID[ix + iy*nbXquads].push_back(i);
	}
}


void pointSet::fillBndQuads()
{
	if (bndQuadID!=NULL)
		delete[] bndQuadID;
	bndQuadID = new vector<int>[nbXquads*nbYquads];
	int ix,iy;
	for (int i=0; i<NbBndPoints; i++)
	{
		getQuad(*pData[i], ix, iy);
		bndQuadID[ix + iy*nbXquads].push_back(i);
	}
}



const vector<int>& pointSet::getQuad(int ix, int iy, bool inner) const
{
	if (inner)
	{
		return quadID[ix + iy*nbXquads];
	}
	else
	{
		return bndQuadID[ix + iy*nbXquads];
	}
}




void pointSet::getQuad(const point& P, int& ix, int& iy) const
{
	ix = max(min(int( (P.X-Bbox[0]) / (slFactorForQuads*SL) ),nbXquads-1),0);
	iy = max(min(int( (P.Y-Bbox[1]) / (slFactorForQuads*SL) ),nbYquads-1),0);
}



void pointSet::findNeighbours(	const point& P, vector<int>& Neighb) const
{
	vector<int> Vn;
	findPossibleNeighbours(P, Vn);
	int NP = Vn.size();
	Neighb.clear();
	double sl2 = pow(getSL(),2);
	for ( int i=0; i<NP; i++)
	{
		double h = P.X - x(Vn[i]);
		double k = P.Y - y(Vn[i]);
		if ( h*h + k*k < sl2 )
			Neighb.push_back(Vn[i]);
	}
}


void pointSet::findPossibleNeighbours(const point& P, vector<int>& Neighb, bool innerQuad) const
{
	Neighb.clear();
	int ix, iy;
	getQuad(P, ix, iy);
	vector<int> V = getQuad(ix,iy,innerQuad);
	for (int i=0; i<V.size(); i++)
	{
		Neighb.push_back(V[i]);
	}
	if (ix<nbXquads-1)
	{
		V = getQuad(ix+1,iy,innerQuad);
		for (int i=0; i<V.size(); i++)
			Neighb.push_back(V[i]);
		if (iy<nbYquads-1)
		{
			V = getQuad(ix+1,iy+1,innerQuad);
			for (int i=0; i<V.size(); i++)
				Neighb.push_back(V[i]);
		}
		if (iy>0)
		{
			V = getQuad(ix+1,iy-1,innerQuad);
			for (int i=0; i<V.size(); i++)
				Neighb.push_back(V[i]);
		}

	}
	if (iy<nbYquads-1)
	{
		V = getQuad(ix,iy+1,innerQuad);
		for (int i=0; i<V.size(); i++)
			Neighb.push_back(V[i]);
	}
	if (ix>0)
	{
		V = getQuad(ix-1,iy,innerQuad);
		for (int i=0; i<V.size(); i++)
			Neighb.push_back(V[i]);
		if (iy<nbYquads-1)
		{
			V = getQuad(ix-1,iy+1,innerQuad);
			for (int i=0; i<V.size(); i++)
				Neighb.push_back(V[i]);
		}
		if (iy>0)
		{
			V = getQuad(ix-1,iy-1,innerQuad);
			for (int i=0; i<V.size(); i++)
				Neighb.push_back(V[i]);
		}
	}
	if (iy>0)
	{
		V = getQuad(ix,iy-1,innerQuad);
		for (int i=0; i<V.size(); i++)
			Neighb.push_back(V[i]);
	}
}

void pointSet::findPossibleNeighbours(int i, vector<int>& Neighb, bool innerQuad) const
{
	findPossibleNeighbours(*pData[i], Neighb, innerQuad);
};


void pointSet::addScalarField(scalarField* sf)
{
	sF.push_back(sf);
}


void pointSet::removeScalarField(scalarField* sf)
{
	vector<scalarField*>::iterator it;
	for (it=sF.begin(); it<sF.end(); it++)
		if (sf == *it)
		{
			sF.erase(it);
			return;
		}
	throw mfpmExcept(10);
}




void pointSet::updateNormals()
{
	vector<int> neighb;
	double nx, ny;
	for (int i=0; i<NbBndPoints; i++)
	{
		if (pData[i]->Bnd <= inactiveFixedWall)
			continue;
		findPossibleNeighbours(i, neighb, false);
		approximatedNormals(*this, neighb, i, nx, ny);
		pData[i]->Nx = nx;
		pData[i]->Ny = ny;
	}

}


void pointSet::adaptSmoothingLength()
{
	int maxSize = 0;
	int minSize = 10000;
	int sum = 0;
	int posMax = -1;
	int posMin = -1;
	for (int i=0; i<NP; i++)
	{
		int n = neighbours[i].size();
		if (minSize > n)
		{
			minSize = n;
			posMin  = i;
		}
		if (maxSize < n)
		{
			maxSize = n;
			posMax  = i;
		}
		sum += n;
	}
	sum /= NP;
	cout << "Size of neighbourhood:  min: " << minSize << "(" << posMin << ") - max: "
	     << maxSize << "(" << posMax << ") - mean: " << sum << endl;
}



int pointSet::updateNeighbours(int minNbP, int maxNbP)
{
	fillQuads();
	if (neighbours != NULL)
		delete[] neighbours;
	neighbours = new vector<int>[NP];
	vector<int> Neighb;
	int maxNb = 0;
	double sl2 = SL*SL;
	double slLoc = sl2;
	int countHighSLchange = 0;
	for (int i=0; i<NP; i++)
	{
		sl2 = pow(getSL(i),2);
		slLoc = sl2;
		// Min: 16, Max:30
		while (neighbours[i].size()<minNbP || neighbours[i].size()>maxNbP)
		{
			//~ if (slLoc > 4*sl2)
			//~ {
			//~ countHighSLchange++;
			//~ break;
			//~ }
			neighbours[i].clear();
			findPossibleNeighbours(*pData[i], Neighb);
			int Vsize = Neighb.size();
			for (int j=0; j<Vsize; j++)
			{
				double dx = pData[Neighb[j]]->X - pData[i]->X;
				double dy = pData[Neighb[j]]->Y - pData[i]->Y;
				if ( dx*dx + dy*dy < slLoc )
				{
					neighbours[i].push_back(Neighb[j]);
				}
			}
			if (neighbours[i].size()<minNbP)
				slLoc *= 1.2;
			else
				slLoc *= 0.8;
		}
		if (slLoc<0.5*sl2)
		{
			countHighSLchange++;
		}
		if (slLoc<0.1*sl2)
		{
			countHighSLchange++;
		}


	}
	return countHighSLchange;
}


int pointSet::updateNeighboursRemoveBad(int minNbP, int maxNbP)
{
	fillQuads();
	if (neighbours != NULL)
		delete[] neighbours;
	neighbours = new vector<int>[NP];
	vector<int> Neighb;
	int maxNb = 0;
	double sl2 = SL*SL;
	double slLoc = sl2;
	int countHighSLchange = 0;
	bool removeBadPoint[NP];
	for (int i=0; i<NP; i++)
	{
		sl2 = pow(getSL(i),2);
		slLoc = sl2;
		removeBadPoint[i] = false;
		// Min: 16, Max:30
		while (neighbours[i].size()<minNbP || neighbours[i].size()>maxNbP)
		{
			if (slLoc > 6*sl2)
			{
				countHighSLchange++;
				removeBadPoint[i] = true;
				break;
			}
			neighbours[i].clear();
			findPossibleNeighbours(*pData[i], Neighb);
			int Vsize = Neighb.size();
			for (int j=0; j<Vsize; j++)
			{
				double dx = pData[Neighb[j]]->X - pData[i]->X;
				double dy = pData[Neighb[j]]->Y - pData[i]->Y;
				if ( dx*dx + dy*dy < slLoc )
				{
					neighbours[i].push_back(Neighb[j]);
				}
			}
			if (neighbours[i].size()<minNbP)
				slLoc *= 1.2;
			else
				slLoc *= 0.8;
		}
		if (slLoc<0.5*sl2)
		{
			countHighSLchange++;
		}
		if (slLoc<0.1*sl2)
		{
			countHighSLchange++;
		}


	}
	removePoints(removeBadPoint);
	return countHighSLchange;
}



void pointSet::findLeftRightBndPoint(int pn, int& Nr, int& Nl) const
{
	// Find a left and right point of a boundary point
	if (pData[pn]->Bnd == 0)
		throw mfpmException(MFPM_EXCEPT,0);
	vector<int> Pn;
	findPossibleNeighbours(pn, Pn, false);
	int size = Pn.size();

	int rightIdx = -1;
	int leftIdx  = -1;
	double rightDist = 1e100;
	double leftDist  = 1e100;

	double tx,nx;
	for (int i=0; i<size; i++)
	{
		if (Pn[i] == pn)
			continue;
		getConeParameters(*pData[pn], Pn[i], tx, nx);
		// TODO: add/remove normal check?
		double dist2 = tx*tx + nx*nx;
		double nDir;
		if (pData[pn]->Bnd == pData[Pn[i]]->Bnd )
			nDir = pData[pn]->Nx*pData[Pn[i]]->Nx + pData[pn]->Ny*pData[Pn[i]]->Ny;
		else
			nDir = 0;
		// TODO: normal check 2
		nDir = 0;
		if (tx>0)
		{
			if (dist2 < rightDist
			        && nDir>=0  )
			{
				rightDist = dist2;
				rightIdx  = Pn[i];
			}
		}
		else
		{
			if (dist2 < leftDist
			        && nDir>=0)
			{
				leftDist = dist2;
				leftIdx  = Pn[i];
			}
		}
	}
	if (rightIdx >= 0)
		Nr = rightIdx;
	else
		Nr = pn;
	if (leftIdx >= 0)
		Nl = leftIdx;
	else
		Nl = pn;

}


void pointSet::removeInnerPoint(int pnb)
{
	if (pnb<NbBndPoints)
		throw mfpmException(MFPM_EXCEPT,1);
	if (pnb>=NP | pnb<0)
		throw mfpmException(MFPM_EXCEPT,2);
	int n = sF.size();
	for (int i=0; i<n; i++)
		sF[i]->removePoint(pnb);
	vector<point*>::iterator it = pData.begin() + pnb;
	delete pData[pnb];
	pData[pnb] = NULL;
	pData.erase(it);
	NP--;
}



void pointSet::removeBndPoint(int pnb)
{
	if (pnb>=NbBndPoints | pnb<0)
		throw mfpmException(MFPM_EXCEPT,2);
	int n = sF.size();
	for (int i=0; i<n; i++)
		sF[i]->removePoint(pnb);
	vector<point*>::iterator it = pData.begin() + pnb;
	delete pData[pnb];
	pData[pnb] = NULL;
	pData.erase(it);
	NP--;
	NbBndPoints--;
}



void pointSet::removePoints(bool toRem[])
{
	int count = 0;
	int n = NP;
	for (int i=0; i<n; i++)
	{
		if (toRem[i])
			if (count<NbBndPoints)
				removeBndPoint(count);
			else
				removeInnerPoint(count);
		else
			count++;
	}
}




void pointSet::addInnerPoint(const point& pnt)
{
	int n = sF.size();
	double val;
	//fillQuads();
	// Add points to scalarFields (without value)
	for (int i=0; i<n; i++)
	{
		sF[i]->addInnerPoint(pnt);
	}
	// Add point to mesh
	point* padd = new point(pnt);
	padd->meshID = meshNumber;
	padd->SL     = getSL(*padd);
	pData.push_back(padd);
	NP++;

}





void pointSet::getConeParameters(const point& P, const point& Pi, double& tau1, double& tau2) const
{
	// calculates parametric representation w.r.t. tangent T and normal N, i.e.
	// P_i = P + tau1*T + tau2*N
	double x = P.X;
	double y = P.Y;
	double nx = P.Nx;
	double ny = P.Ny;
	double xc = Pi.X;
	double yc = Pi.Y;

	tau1 = -ny*(xc-x) + nx*(yc-y);
	tau2 = nx*(xc-x) + ny*(yc-y);
}


void pointSet::getConeParameters(const point& P, int i, double& tau1, double& tau2) const
{
	getConeParameters(P, *pData[i], tau1, tau2);
}

int pointSet::nbPoints() const
{
	return NP;
}


double& pointSet::x(int i)
{
	return pData[i]->X;
}

double& pointSet::y(int i)
{
	return pData[i]->Y;
}

double pointSet::x(int i) const
{
	return pData[i]->X;
}

double pointSet::y(int i) const
{
	return pData[i]->Y;
}

Array pointSet::x() const
{
	Array res(NP);
	for (int i=0; i<NP; i++)
		res(i) = pData[i]->X;
	return res;
}

Array pointSet::y() const
{
	Array res(NP);
	for (int i=0; i<NP; i++)
		res(i) = pData[i]->Y;
	return res;
}


const point& pointSet::P(int i) const
{
	return *pData[i];
}

point& pointSet::P(int i)
{
	return *pData[i];
}

int pointSet::getNIDX(int i) const
{
	return neighbours[i].size();
}



int pointSet::getIDX(int i, int j) const
{
	if (i<0 || i>=NP)
		throw mfpmExcept(10000);
	int n = neighbours[i].size();
	if (j<0 || j>n)
		throw mfpmExcept(10000);
	int res = neighbours[i][j];
	return res;
}


int pointSet::getBND(int i) const
{
	return pData[i]->Bnd;
}

Array pointSet::getBND() const
{
	Array A(NP);
	for(int i=0; i<NP; i++)
		A(i) = getBND(i);
	return A;
}


double pointSet::getSL(const point& P) const
{
	if (SLrefiner == NULL)
		return SL;
	return (*SLrefiner)(P.X, P.Y);

}

double pointSet::getSL(int i) const
{
	if (i>=0 && i<NP)
	{
		return getSL(*pData[i]);
	}
	else
		throw mfpmExcept(2);
}


double pointSet::getSL() const
{
	return SL;
}

double pointSet::getMinSL() const
{
	return SLmin;
}

double pointSet::getMaxSL() const
{
	return SL;
}


int pointSet::getMaxPoints() const
{
	int maxNb = 0;
	for (int i=0; i<NP; i++)
		if (maxNb<neighbours[i].size())
			maxNb = neighbours[i].size();
	return maxNb;
}

double pointSet::getNormal(int comp, int pnt) const
{
	if (comp==0)
		return pData[pnt]->Nx;
	else if (comp == 1)
		return pData[pnt]->Ny;
	else
		throw mfpmException(MFPM_EXCEPT, 3);
}

double& pointSet::getNormal(int comp, int pnt)
{
	if (comp==0)
		return pData[pnt]->Nx;
	else if (comp == 1)
		return pData[pnt]->Ny;
	else
		throw mfpmException(MFPM_EXCEPT, 3);
}

Array pointSet::getNormal(int comp) const
{
	Array res(NP);
	if (comp==0)
		for (int i=0; i<NP; i++)
			res(i) = pData[i]->Nx;
	else if (comp == 1)
		for (int i=0; i<NP; i++)
			res(i) = pData[i]->Ny;
	else
		throw mfpmException(MFPM_EXCEPT, 3);
	return res;

}

void pointSet::setNormal(int direct, int pnt, double val)
{
	if ( pnt >= NP || pnt < 0 )
		throw mfpmExcept(20);
	if (direct==0)
		pData[pnt]->Nx = val;
	else if (direct == 1)
		pData[pnt]->Ny = val;
	else
		throw mfpmException(MFPM_EXCEPT, 3);
}


void pointSet::setBND(int i, int val)
{
	pData[i]->Bnd = val;
}



int pointSet::getNbBndPoints() const
{
	return NbBndPoints;
}


double pointSet::dist(int p1, int p2) const
{
	return pData[p1]->dist(*pData[p2]);
}


void pointSet::operator=(const pointSet& ps)
{
	if (&ps == this)
		throw mfpmExcept(24);
	int n = sF.size();
	for (int i=0; i<n; i++)
	{
		sF[i]->changeMesh(ps);
	}
	copyPointSet(ps);
}

void pointSet::copyMeshWithUpdate(const pointSet& ps)
{
	if (&ps == this)
		throw mfpmExcept(24);
	int n = sF.size();
	for (int i=0; i<n; i++)
	{
		sF[i]->changeMesh(ps);
	}
	copyPointSet(ps);
}

void pointSet::copyMesh(const pointSet& ps)
{
	if (&ps == this)
		throw mfpmExcept(24);
	copyPointSet(ps);
}

void pointSet::movePoints(vectorField& V)
{
	for (int i=0; i<NP; i++)
	{
		pData[i]->X += V.Xcomp(i);
		pData[i]->Y += V.Ycomp(i);
	}
}

void pointSet::movePoints(const Array& Vx, const Array& Vy)
{
	for (int i=0; i<NP; i++)
	{
		pData[i]->X += Vx(i);
		pData[i]->Y += Vy(i);
	}
}



void pointSet::clearRefiner()
{
	SLrefiner = NULL;
}


void pointSet::setRefiner(refinement& SF)
{
	SLrefiner = &SF;
	for(int i=0; i<NP; i++)
	{
		pData[i]->SL = SF(pData[i]->X,pData[i]->Y);
	}
}



void pointSet::printParameters()
{
	PS_parameter.print();
}

void pointSet::printGeoInfo()
{
	ofstream fil(name.c_str());
	fil << geometry.nbLines << endl;
	fil << geometry.nbDataPoints << endl;
	fil << geometry.nbPoints << endl;
	for (int i=0; i<geometry.nbLines; i++)
	{
		fil << geometry.active[i];
		for (int j=0; j<3; j++)
			fil<<geometry.line[i][j];
		fil << geometry.normalX[i];
		fil << geometry.normalY[i];
	}

	for (int i=0; i<geometry.nbPoints; i++)
	{
		fil << geometry.P[i].X;
		fil << geometry.P[i].Y;
		fil << geometry.P[i].Bnd;
		fil << geometry.P[i].id;
		fil << geometry.P[i].Nx;
		fil << geometry.P[i].Ny;
		fil << geometry.P[i].active;

	}

	for (int i=0; i<geometry.nbDataPoints; i++)
	{
		fil << geometry.dataPoints[i].X;
		fil << geometry.dataPoints[i].Y;
		fil << geometry.dataPoints[i].Bnd;
		fil << geometry.dataPoints[i].id;
		fil << geometry.dataPoints[i].Nx;
		fil << geometry.dataPoints[i].Ny;
		fil << geometry.dataPoints[i].active;
	}

}


void pointSet::printPointData()
{
	ofstream fil(name.c_str());
	fil << NP;
	fil << NbBndPoints;
	fil << SL;
	fil << nbXquads;
	fil << nbYquads;

	for (int i=0; i<NP; i++)
	{
		fil << pData[i]->X;
		fil << pData[i]->Y;
		fil << pData[i]->id;
		fil << pData[i]->Bnd;
		fil << pData[i]->Nx;
		fil << pData[i]->Ny;
		fil << pData[i]->active;
	}

	for (int i=0; i<nbXquads*nbYquads; i++)
	{
		for (int j=0; j<quadID[i].size(); j++)
		{
			fil << quadID[i][j];
		}
		for (int j=0; j<bndQuadID[i].size(); j++)
		{
			fil << bndQuadID[i][j];
		}
	}

	for (int i=0; i<NP; i++)
	{
		for (int j=0; j<neighbours[i].size(); j++)
		{
			fil << neighbours[i][j];
		}
	}


}


int pointSet::differentBoundaries()
{
	int bnds = 1;
	int bID[500];
	bID[0] = pData[0]->Bnd;
	for (int i=0; i<getNbBndPoints(); i++)
	{
		bool found = false;
		for (int j=0; j<bnds; j++)
		{
			if ( pData[i]->Bnd == bID[j] )
			{
				found=true;
				break;
			}
		}
		if ( !found )
		{
			bID[bnds] = pData[i]->Bnd;
			bnds++;
		}
	}
	return bnds;
}



void pointSet::write_ADOLC_Array(
    int bnd_block_size,			// size of one BND block
    int inner_block_size,		// size of block for inner points
    int nb_bnd_blocks,			// Nb. of boundary blocks, i.e. nb of different boundaries
    int max_nb_neighbours,		// Max. nb. of neighbours (e.g. 40?)
    int* bnd_block_list,		// Boundary list, i.e. for boundaries 10,20,30 bnd_block_list[0]=10, bnd_block_list[0]=20...
    double* X_arr,				// Array with x-coord
    double* Y_arr,				// Array with y-coord
    int*    Plist				// Pointer list of neighbours organised as Plist=[ |BND1| |BND2| ... |INNER| ]
)								// |BND1| = [ Neighb(x1) | Neighb(x2) | .... | -1 -1 ... ]
// Neighb(x1) = all neighbours of point x1, diagonal element at front
{
	// Clear X_arr, Y_arr
	for (int i=0; i<bnd_block_size*nb_bnd_blocks+inner_block_size; i++)
	{
		X_arr[i] = 0.0;
		Y_arr[i] = 0.0;
	}

	// Fill boundary points
	for( int i=0; i<NbBndPoints; i++ )
	{
		X_arr[i] = pData[i]->X;
		Y_arr[i] = pData[i]->Y;
	}

	// Fill inner points
	for( int i=NbBndPoints; i<NP; i++ )
	{
		X_arr[i] = pData[i]->X;
		Y_arr[i] = pData[i]->Y;
	}

	// Clear Plist
	for(int i=0; i<NP*max_nb_neighbours; i++)
		Plist[i] = -1;

	// Generate point list array for each BND block
	for( int bnd_block=0; bnd_block<nb_bnd_blocks; bnd_block++)
	{
		int counter=bnd_block*bnd_block_size*max_nb_neighbours;
		for (int i=0; i<NbBndPoints; i++)
		{
			if( pData[i]->Bnd == bnd_block_list[bnd_block] )
			{
				for (int j=0; j<getNIDX(i); j++)
				{
					if( getIDX(i,j) == i )
					{
						Plist[counter] = getIDX(i,j);
						counter++;
						break;
					}
				}
				for (int j=0; j<getNIDX(i); j++)
				{
					if ( getIDX(i,j) != i)
					{
						Plist[counter] = getIDX(i,j);
						counter++;
					}
				}
				for (int j=getNIDX(i); j<max_nb_neighbours; j++)
				{
					Plist[counter] = -1;
					counter++;
				}
				cout << "P" << i << " : " << counter << endl;

			}
		}
	}

	// Generate point list array for inner points
	int counter = nb_bnd_blocks*bnd_block_size*max_nb_neighbours;
	for( int i=NbBndPoints; i<NP; i++)
	{
		for (int j=0; j<getNIDX(i); j++)
		{
			if( getIDX(i,j) == i )
			{
				Plist[counter] = getIDX(i,j);
				counter++;
				break;
			}
		}
		for (int j=0; j<getNIDX(i); j++)
		{
			if ( getIDX(i,j) != i)
			{
				Plist[counter] = getIDX(i,j);
				counter++;
			}
		}
		for (int j=getNIDX(i); j<max_nb_neighbours; j++)
		{
			Plist[counter] = -1;
			counter++;
		}
		cout << "P" << i << " : " << counter << endl;
	}
}

