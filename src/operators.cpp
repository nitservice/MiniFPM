#include "operators.h"
#include <omp.h>
#include <cmath>
#include "Array.h"
#include "pointSet.h"
#include "scalarField.h"
#include "vectorField.h"
#include "primitive_operators.h"


bool smoothFieldWarningPrinted = false;

mfpmMatrix laplace(const pointSet& mesh, bool derivBnd)
{
	mfpmMatrix mat(mesh);
	int nbNotCon;
	if (derivBnd)
		nbNotCon = 0;
	else
		nbNotCon = mesh.getNbBndPoints();

	Array tmp(200);
	int i, info;

//#pragma omp parallel for default(shared) private(i, tmp, info)
	for (int i=0; i<nbNotCon; i++)
	{
		tmp.setSize(mesh.getNIDX(i));
		mat.setRow(i, tmp);
	}

//#pragma omp parallel for default(shared) private(i, tmp, info)
	for (i=nbNotCon; i<mesh.nbPoints(); i++)
	{
		//Array tmp(mesh.getNIDX(i));
		tmp.setSize(mesh.getNIDX(i));
		info = laplaceMatEntry(mesh, i, tmp);
		if (info != 0)
			cout << "!" << info << "!" <<  endl;
		mat.setRow(i, tmp);
	}


	return mat;

}


mfpmMatrix laplace4(const pointSet& mesh, bool derivBnd)
{
	mfpmMatrix mat(mesh);
	int nbNotCon;
	if (derivBnd)
		nbNotCon = 0;
	else
		nbNotCon = mesh.getNbBndPoints();

//#pragma omp parallel for default(shared)
	for (int i=0; i<nbNotCon; i++)
	{
		Array tmp(mesh.getNIDX(i));
		mat.setRow(i, tmp);
	}
	int errCnt = 0;

//#pragma omp parallel for default(shared)
	for (int i=nbNotCon; i<mesh.nbPoints(); i++)
	{
		Array tmp(mesh.getNIDX(i));
		int info = laplaceOrder4MatEntry(mesh, i, tmp);
		if (info != 0)
			errCnt++;
		mat.setRow(i, tmp);
	}
	//cout << "Errors: " << errCnt << " / " << mesh.nbPoints()-nbNotCon << endl;
	return mat;

}


mfpmMatrix biLaplace(const pointSet& mesh)
{
	mfpmMatrix mat(mesh);
	for (int i=0; i<mesh.getNbBndPoints(); i++)
	{
		Array tmp(mesh.getNIDX(i));
		mat.setRow(i, tmp);
	}
	for (int i=mesh.getNbBndPoints(); i<mesh.nbPoints(); i++)
	{
		Array tmp(mesh.getNIDX(i));
		int info = biLaplaceMatEntry(mesh, i, tmp);
		if (info != 0)
			cout << "!" << info << "!" <<  endl;
		mat.setRow(i, tmp);
	}
	return mat;

}


mfpmMatrix gradNormal(const pointSet& mesh)
{
	mfpmMatrix mat(mesh);

	for (int i=0; i<mesh.nbPoints(); i++)
	{
		if (mesh.getBND(i) != 0)
		{
			Array tmp(mesh.getNIDX(i));
			neumannMatEntry(mesh, i, tmp);
			mat.setRow(i, tmp);
		}
	}
	return mat;
}


void gradient(const pointSet& mesh, mfpmMatrix& Gx, mfpmMatrix& Gy, int order)
{
//	#pragma omp parallel for default(shared)
	for (int i=0; i<mesh.nbPoints(); i++)
	{
		Array tmpX(mesh.getNIDX(i));
		Array tmpY(mesh.getNIDX(i));
		int locOrder = order;
		int info = gradMatEntry(mesh, i, tmpX, tmpY, order);
		while (info != 0)
		{
			locOrder--;
			info = gradMatEntry(mesh, i, tmpX, tmpY, locOrder);
			if (info == 0)
				cout << "Gradient expects order " << order << " but reduces to order " << locOrder << endl;
		}
		Gx.setRow(i, tmpX);
		Gy.setRow(i, tmpY);
	}
}


mfpmMatrix ddt(const pointSet& mesh, double tau, scalarField& oldTime)
{
	if (&mesh != oldTime.getMesh())
		throw mfpmExcept(25);
	mfpmMatrix mat(mesh);
	double** data = mat.getData();
	Array* rhs = mat.getRhs();
	for (int i=0; i<mesh.nbPoints(); i++)
	{
		for (int j=0; j<mesh.getNIDX(i); j++)
			if (mesh.getIDX(i,j) == i)
				data[i][j] = 1.0 / tau;
			else
				data[i][j] = 0.0;
		double oldVal = oldTime(i);
		rhs->operator()(i) = oldVal / tau;
	}
	return mat;
}


mfpmMatrix isource(const pointSet& mesh)
{
	mfpmMatrix mat(mesh);
	double** data = mat.getData();
	Array* rhs = mat.getRhs();
	for (int i=0; i<mesh.nbPoints(); i++)
	{
		for (int j=0; j<mesh.getNIDX(i); j++)
			if (mesh.getIDX(i,j) == i)
				data[i][j] = 1.0;
			else
				data[i][j] = 0.0;
	}
	return mat;
}



void explDerivatives(const scalarField& f, vectorField* grad, scalarField* lapF, scalarField* biLapF, bool inner)
{
	const pointSet& mesh = *f.getMesh();
	double a[15];
	for (int i=0; i<mesh.nbPoints(); i++)
	{
		int order = interpolationConstants(f, i, a, inner);
		double x = mesh.P(i).X;
		double y = mesh.P(i).Y;
		if (lapF != NULL)
		{
			if (order < 2)
				throw mfpmExcept(26);
			(*lapF)(i) = 2*a[3] + 6*a[6]*x + 2*a[7]*y + 12*a[10]*x*x + 6*a[11]*x*y + 2*a[12]*y*y
			             + 2*a[5] + 2*a[8]*x + 6*a[9]*y + 2*a[12]*x*x + 6*a[13]*x*y + 12*a[14]*y*y;
		}

		if (biLapF != NULL)
		{
			if (order < 4)
				throw mfpmExcept(26);
			(*biLapF)(i) = 24*a[10] + 8*a[12] + 24*a[14];
		}

		if (grad != NULL)
		{
			if (order < 1)
				throw mfpmExcept(26);
			grad->Xcomp()(i) = a[1] + 2*a[3]*x + a[4]*y + 3*a[6]*x*x + 2*a[7]*x*y + a[8]*y*y + 4*a[10]*x*x*x
			                   + 3*a[11]*x*x*y + 2*a[12]*x*y*y + a[13]*y*y*y;

			grad->Ycomp()(i) = 	a[2] + a[4]*x + 2*a[5]*y + a[7]*x*x + 2*a[8]*x*y + 3*a[9]*y*y + a[11]*x*x*x
			                    + 2*a[12]*x*x*y + 3*a[13]*x*y*y + 4*a[14]*y*y*y;
		}
	}
}


scalarField smoothField(scalarField& sf, double alpha, double smoothingLength)
{
	pointSet& mesh = *sf.getMesh();
	double h;
	if (smoothingLength <= 0)
		h = mesh.getMinSL();
	else
		h = smoothingLength;
	scalarField tmp(mesh);
	double dist;

	h = h*alpha;

	if (h>mesh.getMaxSL())
	{
		if (!smoothFieldWarningPrinted)
		{
			cout << "!!!" << endl;
			cout << "Warning: Smoothing factor for '" << sf.name << "' too large -> h is reduced to hmax!" << endl;
			cout << "!!!" << endl;
			smoothFieldWarningPrinted = true;
		}
		//h = mesh.getMaxSL();
	}

	for (int i=0; i<mesh.nbPoints(); i++)
	{
		vector<int> neighb;
		mesh.findPossibleNeighbours(i, neighb);
		double weight;
		double sum = 0.0;
		for (int j=0; j<neighb.size(); j++)
		{
			// Perhaps also "dist2" gives good results?
			dist = mesh.P(i).dist(mesh.P(neighb[j])) / (h);
			if ( dist<1 )
				weight = 1.0 - 6.0*dist*dist + 8.0*dist*dist*dist - 3.0*dist*dist*dist*dist;
			else
				weight = 0.0;
			sum += weight;
			tmp(i) += weight * sf(neighb[j]);
		}
		tmp(i) = tmp(i) / sum;

	}

	return tmp;


}




scalarField smoothField_full(scalarField& sf, double alpha, double smoothingLength)
{
	pointSet& mesh = *sf.getMesh();
	double h;
	if (smoothingLength <= 0)
		h = mesh.getMinSL();
	else
		h = smoothingLength;
	scalarField tmp(mesh);
	double dist;

	h = h*alpha;

	int n = mesh.nbPoints();
	for (int i=0; i<n; i++)
	{
		double weight;
		double sum = 0.0;
		tmp(i) = 0.0;
		for (int j=0; j<n; j++)
		{
			// Perhaps also "dist2" gives good results?
			dist = mesh.P(i).dist(mesh.P(j)) / (h);
			if ( dist<1 )
				weight = 1.0 - 6.0*dist*dist + 8.0*dist*dist*dist - 3.0*dist*dist*dist*dist;
			else
				weight = 0.0;
			sum += weight;
			tmp(i) += weight * sf(j);
		}
		if (sum==0)
			throw mfpmExcept(100000);
		tmp(i) = tmp(i) / sum;

	}

	return tmp;


}



mfpmMatrix divEtaGrad(const scalarField& mu)
{
	const pointSet& mesh = *mu.getMesh();
	int NP = mesh.nbPoints();
	bool interfacialPoint[NP];

	mfpmMatrix divMuGrad(mesh);

	/********************************
	* Determine interfacial points
	********************************/

	double SL = 0.4*mesh.getSL();
	double SL2 = SL*SL;
	int nbInterPnt = 0;
	for (int i=0; i<mesh.nbPoints(); i++)
	{
		double muMax = -1e100;
		double muMin = 1e100;
		for (int j=0; j<mesh.getNIDX(i); j++)
		{
			if (mu(mesh.getIDX(i,j)) > muMax)
				muMax = mu(mesh.getIDX(i,j));
			if (mu(mesh.getIDX(i,j)) < muMin)
				muMin = mu(mesh.getIDX(i,j));
		}
		double avgMu = (muMax+muMin)/2.0;

		double xp = mesh.x(i);
		double yp = mesh.y(i);

		interfacialPoint[i] = false;

		//if (muMin > 0.1*muMax ) continue;

		for (int j=0; j<mesh.getNIDX(i); j++)
		{
			double dist = pow(mesh.x(mesh.getIDX(i,j))-xp,2) + pow(mesh.y(mesh.getIDX(i,j))-yp,2);
			// FIXME dist is set to fixed value!
			dist = 0.0;
			if ( mu(i) >= avgMu && mu(mesh.getIDX(i,j)) <= avgMu &&  dist < SL2 )
			{
				interfacialPoint[i] = true;
				nbInterPnt++;
				break;
			}
			else if ( mu(i) < avgMu && mu(mesh.getIDX(i,j)) >= avgMu &&  dist < SL2 )
			{
				interfacialPoint[i] = true;
				nbInterPnt++;
				break;
			}
		}
	}
	//cout << "Number of interfacial points:  " << nbInterPnt  <<  endl;

	/********************************
	* Build laplace matrix for left and right interface area
	********************************/
	for (int i=mesh.getNbBndPoints(); i<mesh.nbPoints(); i++)
	{
		if ( !interfacialPoint[i] )
		{
			//~ Array active(mesh.getNIDX(i));
			//~ for (int j=0; j<mesh.getNIDX(i); j++)
			//~ {
			//~ // TODO: adapt for arbitrary mu's!!!
			//~ if ( mu(mesh.getIDX(i,j)) > 1.0 )
			//~ active(j) = 1.0;
			//~ else
			//~ active(j) = 0.0;
			//~ }
			//~ if (mu(i)<1.0)
			//~ active = 1.0-active;

			Array tmp(mesh.getNIDX(i));
			int info = laplaceMatEntry(mesh, i, tmp);
			if (info != 0)  cout << "!" << info << "!" <<  endl;
			divMuGrad.setRow(i, mu(i)*tmp);
		}
	}

	/*********************************
	 * Approximate Laplacian in interfacial region
	 ********************************/
	for (int part=mesh.getNbBndPoints(); part<mesh.nbPoints(); part++)
	{
		if (!interfacialPoint[part] )
			continue;

		int info;

		double sl = mesh.getSL();
		double alpha = 0.5;
		double dx = alpha*sl;
		double dy = alpha*sl;

		double xp = mesh.x(part);
		double yp = mesh.y(part);
		double xn[4], yn[4], lright, lleft;
		// Define approximaion points
		xn[0] = xp + dx;
		xn[1] = xp - dx;
		xn[2] = xp;
		xn[3] = xp;

		yn[0] = yp;
		yn[1] = yp;
		yn[2] = yp + dy;
		yn[3] = yp - dy;


		Array entLeft(mesh.getNIDX(part));
		Array entRight(mesh.getNIDX(part));
		Array entDiag(mesh.getNIDX(part));

		// Derive second derivative in x-direction
		info = approximationXYMatEntry(mesh, part, xn[0], yn[0], entRight,1);
		info = approximationXYMatEntry(mesh, part, xn[1], yn[1], entLeft,1);

		// Derive second derivative in x-direction with order 2
		//~ Array entRight2(mesh.getNIDX(part));
		//~ Array entLeft2(mesh.getNIDX(part));
		//~ info = approximationXYMatEntry(mesh, part, xn[0], yn[0], entRight2,2);
		//~ info = approximationXYMatEntry(mesh, part, xn[1], yn[1], entLeft2,2);

		// Approximate lambda plus und minus and diagonal element
		lright  = 0;
		lleft = 0;
		for (int i=0; i<mesh.getNIDX(part); i++)
		{
			lright += entRight(i)*mu(mesh.getIDX(part,i));
			lleft += entLeft(i)*mu(mesh.getIDX(part,i));
			if (mesh.getIDX(part,i) == part)
				entDiag(i) = 1.0;
			else
				entDiag(i) = 0.0;
		}

		divMuGrad.setRow(part, (lright*entRight - (lright+lleft)*entDiag + lleft*entLeft)/(dx*dx) );



		// Derive second derivative in y-direction
		info = approximationXYMatEntry(mesh, part, xn[2], yn[2], entRight,1);
		info = approximationXYMatEntry(mesh, part, xn[3], yn[3], entLeft,1);

		// Derive second derivative in x-direction with order 2
		//~ info = approximationXYMatEntry(mesh, part, xn[2], yn[2], entRight2,2);
		//~ info = approximationXYMatEntry(mesh, part, xn[3], yn[3], entLeft2,2);

		// Approximate lambda plus und minus and diagonal element
		lright  = 0;
		lleft = 0;
		for (int i=0; i<mesh.getNIDX(part); i++)
		{
			lright += entRight(i)*mu(mesh.getIDX(part,i));
			lleft += entLeft(i)*mu(mesh.getIDX(part,i));
			if (mesh.getIDX(part,i) == part)
				entDiag(i) = 1.0;
			else
				entDiag(i) = 0.0;
		}

		divMuGrad.addRow(part, (lright*entRight - (lright+lleft)*entDiag + lleft*entLeft)/(dy*dy) );
	}
	return divMuGrad;
}

