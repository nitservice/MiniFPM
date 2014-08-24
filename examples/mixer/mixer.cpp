#include <cmath>
//#include <iostream>
#include <cstdio>

#include "pointSet.h"
#include "scalarField.h"
#include "postResult.h"
#include "primitive_operators.h"
#include "Array.h"
#include "mfpmMatrix.h"
#include "operators.h"

#include "matrixUmfpack.h"

using namespace std;

int main(int argc, char** argv)
{
	char ct[100];
	mfpmParameters para;
	para.print();

	pointSet mesh("geometry.msh", -1, 200);
	postResult post(mesh);


	post.writeGeometryVTK("Results/res_00000.vtk");

	int n = mesh.nbPoints();

	double tau = 0.015;
	int    nt  = 150;
	double mu  = 1e-1;

	try
	{

		scalarField u(mesh, "U");
		scalarField v(mesh, "V");
		scalarField p(mesh, "P");
		scalarField pcorr(mesh, "Pcorr");
		scalarField divNew(mesh, "divNew");
		scalarField divOld(mesh, "divOld");
		vectorField velo(mesh, "velo");
		scalarField& velox = velo.Xcomp();
		scalarField& veloy = velo.Ycomp();
		vectorField normals(mesh, "Normals");
		scalarField KOB(mesh, "BndID");
		scalarField pHyd(mesh, "pHyd");
		scalarField partID(mesh, "partID");
		scalarField vorticity(mesh, "vorticity");
		post.add(vorticity);
		post.add(partID);
		post.add(p);
		post.add(pHyd);
		post.add(divNew);
		post.add(divOld);
		post.add(velo);
		post.add(normals);
		post.add(KOB);

		u = 0.0;
		v = 0.0;
		p = 0.0;


		for (int i=0; i<mesh.nbPoints(); i++)
		{
			partID(i) = i;
		}


		for (int tl=1; tl<nt; tl++)
		{
			cout << "Time-step " << tl << endl;

			for (int k=0; k<mesh.nbPoints(); k++)
			{
				if (mesh.getBND(k) >= 200 || mesh.getBND(k) == 0)
				{
					if (abs(u(k)*tau) > 2*mesh.getSL(0) || abs(v(k)*tau) > 2*mesh.getSL(0) )
					{
						cout << "Deformation too high -> exit" << endl;
						//throw mfpmExcept(1000);
					}
					mesh.x(k) = mesh.x(k) + u(k)*tau;
					mesh.y(k) = mesh.y(k) + v(k)*tau;
				}
			}

			cout << " -> adapt boundary points" << endl;
			mesh.adaptWallBoundaryPoints();
			mesh.adaptFreeSurfaceBoundaryPoints();

			cout << " -> add inner points" << endl;
			mesh.adaptRemoveInnerPoints();
			int nbAdded;
			nbAdded = mesh.adaptAddInnerPoints();
			cout << "   --> Nb points added: " << nbAdded << endl;

			cout << " -> updateNeighbours" << endl;
			mesh.updateNeighbours();
			mesh.adaptSmoothingLength();

			n = mesh.nbPoints();

			mfpmMatrix A(mesh);
			mfpmMatrix Gx(mesh);
			mfpmMatrix Gy(mesh);
			mfpmMatrix lap = laplace(mesh);

			gradient(mesh,Gx,Gy);


			// Solve for hydrostatic pressure
			cout << " -> build gradient matrices... " << endl;
			Array hpbnd(n);
			double fx = 0.0;
			double fy = -10.0;
			for (int i=0; i<n; i++)
				hpbnd(i) = (fx*mesh.getNormal(0,i) + (fy)*mesh.getNormal(1,i));
			A = lap;
			A.setNeumann(1, hpbnd);
			A.setNeumann(100, hpbnd);
			A.setDirichlet(200, 0.0);
			A.setNeumann(201, hpbnd);
			A.solveUMFPACK(pHyd);



			A = ddt(mesh, tau, u) - mu*lap == fx - Gx*pHyd;
			A.setDirichlet(1, 0.0);
			A.setDirichlet(100, 0.0);
			A.setNeumann(200, 0.0);
			A.setDirichlet(201, sin(tl*tau*2*M_PI));
			A.solveUMFPACK(u);

			A = ddt(mesh, tau, v) - mu*lap == fy - Gy*pHyd;
			A.setDirichlet(1, 0.0);
			A.setDirichlet(100, 0.0);
			A.setNeumann(200, 0.0);
			A.setDirichlet(201, cos(tl*tau*2*M_PI));
			A.solveUMFPACK(v);


			Array divU, gradPx, gradPy;

			divU = Gx*u + Gy*v;

			A = tau*lap == divU;
			A.setNeumann(1, 0.0);
			A.setNeumann(100, 0.0);
			A.setDirichlet(200, 0.0);
			A.setNeumann(201, 0.0);

			A.solveUMFPACK(pcorr);

			gradPx = Gx*pcorr;
			gradPy = Gy*pcorr;

			for (int i=0; i<n; i++)
				if (mesh.getBND(i) == 200
				        || mesh.getBND(i) == 0)
				{
					u(i) = u(i) - tau*gradPx(i);
					v(i) = v(i) - tau*gradPy(i);
				}

			p() = pcorr();

			// Post-processing
			for (int i=0; i<n; i++)
			{
				KOB(i) = mesh.getBND(i);
			}
			velox = u;
			veloy = v;
			vorticity = Gx*v - Gy*u;
			normals.Xcomp() = mesh.getNormal(0);
			normals.Ycomp() = mesh.getNormal(1);
			sprintf(ct, "Results/res_%5.5i.vtk", tl);
			std::string fn(ct);
			post.writeVTK(fn);
		}


		return 0;

	}
	catch (const mfpmException& mex)
	{
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "!! ERROR OCCURRED!" << endl;
		cout << "!! ==============" << endl;
		cout << "!! " << mex.what() << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		post.writeGeometryVTK("test1.vtk");

	}

}



