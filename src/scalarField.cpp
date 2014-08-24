#include "scalarField.h"
#include <fstream>
#include <cmath>
#include "pointSet.h"
#include "primitive_operators.h"

scalarField::scalarField(pointSet &mesh)
{
	this->mesh = &mesh;
	this->name = "tmpField";
	this->data.setSize(mesh.nbPoints());
	mesh.addScalarField(this);
}

scalarField::scalarField(pointSet &mesh, const string &name)
{
	this->mesh = &mesh;
	this->name = name;
	this->data.setSize(mesh.nbPoints());
	mesh.addScalarField(this);
}


scalarField::scalarField(const scalarField& sf)
{
	this->mesh = sf.getMesh();
	this->name = sf.name;
	data = sf();
	mesh->addScalarField(this);
}


scalarField::scalarField()
{
	this->mesh = NULL;
	this->name = "tmpField";
}


scalarField::~scalarField()
{
	if (mesh != NULL)
	{
		mesh->removeScalarField(this);
		mesh = NULL;
		this->name = "";
	}
}


double& scalarField::operator()(int i)
{
	return data(i);
}

double scalarField::operator()(int i) const
{
	return data(i);
}


double scalarField::operator()(double x, double y) const
{
	point pTmp;
	pTmp.X = x;
	pTmp.Y = y;
	return interpolateLS(pTmp);
}


void scalarField::operator=(const scalarField& sf)
{
	pointSet* ps = sf.getMesh();

	if (data.empty())
	{
		this->mesh = sf.getMesh();
		this->name = sf.name;
		this->data = sf();
		mesh->addScalarField(this);
	}
	else
	{
		if (ps == this->mesh)
			data = sf();
		else
		{
			int n = data.getSize();
			for (int i=0; i<n; i++)
			{
				if ( i < mesh->getNbBndPoints() )
					data(i) = sf.interpolate(mesh->P(i));
				else
					data(i) = sf.interpolateLS(mesh->P(i));
			}
		}

	}
}

void scalarField::operator=(const Array& v)
{
	if (data.empty())
	{
		throw "Error scalarField:: operator=: object not initialised!";
	}
	else
	{
		if (v.getSize() == this->getSize())
			data = v;
		else
			throw "Error scalarField:: operator=: lshape != rshape.";
	}
}


void scalarField::operator=(const double val)
{
	if (data.empty())
	{
		throw "Error scalarField:: operator=: object not initialised!";
	}
	else
	{
		data = val;
	}
}





void scalarField::print() const
{
	cout << "ScalarField " << name << ":" << endl;
	for (int i=0; i<getSize(); i++)
		cout << data(i) << endl;
}




void scalarField::removePoint(int i)
{
	data.remove(i);
}


void scalarField::addPoint(int i, double val)
{
	data.add(i, val);
}


void scalarField::addPoint(int i)
{
	data.add(i);
}

void scalarField::addInnerPoint(const point& P)
{
	double val = interpolateLS(P);
	data.append(val);
}


void scalarField::addBndPoint(int pos, const point& P)
{
	double val = interpolate(P, false);
	data.add(pos);
	data(pos) = val;
}



double scalarField::interpolate(const point& P, bool inner) const
{
	vector<int> neighb;
	mesh->findPossibleNeighbours(P, neighb, false);
	int n = neighb.size();
	if (n == 0)
	{
		// TODO: Validate zero interpolation!
		return 0.0;
	}
	double w[n];
	double sum = 0;
	double fac;
	double h = mesh->getSL(P);
	fac=1.0/h*h;
	// Build preliminary weights
	for (int i=0; i<n; i++)
	{
		if (P.Bnd == mesh->P(neighb[i]).Bnd)
		{
			double dx = P.dist2(mesh->P(neighb[i])) ;
			// TODO: Adapt exp factor
			w[i] = exp(-300.0*dx * fac );
			sum += w[i];
		}
		else
			w[i] = 0.0;
	}

	if (sum==0)
	{
		if (n>0)
			return data(neighb[0]);
		else
		{
			if (inner)
				throw mfpmExcept(20);
			else
				throw mfpmExcept(21);
		}
	}
	// Normalise (Shepard)
	for (int i=0; i<n; i++)
		w[i] /= sum;
	// Calc. interpolated value
	sum = 0;
	for (int i=0; i<n; i++)
		sum += w[i]*data(neighb[i]);

	return sum;

}



double scalarField::interpolateLS(const point& P, bool forceLowOrder) const
{
	vector<int> neighb;
	mesh->findPossibleNeighbours(P, neighb);
	int m = neighb.size();
	double rhs[m];
	double A[6*m];
	int mi;
	double weight;
	double h2 = pow((mesh->getMaxSL()+mesh->getMinSL())/2,2);
	//double h2 = pow( mesh->getSL(P) , 2);
	int orderFac;
	if (m<3)
	{
		// TODO: Validate zero return!
		return 0.0;
		throw mfpmExcept(20);
	}

	// Check whether interpolation is needed
	double dmin = 1e100;
	double idx_min;
	const double tolerance = 0.05* mesh->getMinSL();
	for (int i=0; i<m; i++)
	{
		double d = P.dist(mesh->P(neighb[i]));
		if ( d < dmin )
		{
			dmin = d;
			idx_min = i;
		}
	}
	if ( dmin < tolerance )
	{
		return value(neighb[idx_min ]);
	}

	if (forceLowOrder | m<6)
		orderFac = 3;
	else
		orderFac = 6;

	for (int i=0; i<m; i++)
	{
		const point& Ptmp = mesh->P(neighb[i]);
		mi = orderFac*i;
		weight = exp(-3.0*P.dist2(Ptmp) / h2);
		A[mi+0] = weight*1.0;
		A[mi+1] = weight*Ptmp.X;
		A[mi+2] = weight*Ptmp.Y;
		if (orderFac>3)
		{
			A[mi+3] = weight*Ptmp.X*Ptmp.Y;
			A[mi+4] = weight*Ptmp.X*Ptmp.X;
			A[mi+5] = weight*Ptmp.Y*Ptmp.Y;
		}
		rhs[i]  = weight*value(neighb[i]);
	}

	int n = orderFac;
	double work[m+n];
	int info, workl = m+n;
	int one = 1;
	char TN[] = "T";
	dgels_(TN, &n, &m, &one, A, &n, rhs, &m, work, &workl, &info);

	if (orderFac>3)
	{
		return rhs[0] + rhs[1]*P.X + rhs[2]*P.Y + rhs[3]*P.X*P.Y
		       + rhs[4]*P.X*P.X + rhs[5]*P.Y*P.Y;
	}
	else
	{
		return rhs[0] + rhs[1]*P.X + rhs[2]*P.Y;
	}




}



void scalarField::changeMesh(const pointSet& meshNew)
{
	int n = meshNew.nbPoints();
	Array tmp(n);
	for (int i=0; i<n; i++)
	{
		if ( i < meshNew.getNbBndPoints() )
		{
			tmp(i) = interpolate(meshNew.P(i), false);
		}
		else
		{
			tmp(i) = interpolateLS(meshNew.P(i));
		}
	}
	data.setSize(n);
	data = tmp;
}



void scalarField::assignBndValue(const scalarField& sf, int bnd)
{
	//mesh->fillBndQuads();
	//sf.getMesh()->fillBndQuads();
	for (int i=0; i<mesh->getNbBndPoints(); i++)
	{
		if (mesh->getBND(i) == bnd)
		{
			data(i) = sf.interpolate(mesh->P(i), false);
		}
	}
}




void scalarField::printData()
{
	ofstream fil(name.c_str());
	for (int i=0; i<data.getSize(); i++)
	{
		fil << data(i) << endl;
	}
}
