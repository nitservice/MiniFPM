#include "primitive_operators.h"
#include <cmath>
#include "pointSet.h"
#include "scalarField.h"
#include "Array.h"

using namespace std;

int laplaceMatEntry(const pointSet& mesh, int pnb, Array& entries)
{
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[7*NP];
	double weight[NP];

	if (NP < 7)
		throw mfpmExcept(20);

	b[0] = -10.0 / pow(mesh.getSL(pnb),2);
	b[1] = 0.0;
	b[2] = 0.0;
	b[3] = 0.0;
	b[4] = 2.0;
	b[5] = 2.0;
	b[6] = 0.0;

	int diagEntry;
	int idx, mi;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = 7*i;
		if (idx!=pnb)
		{
			h = mesh.x(idx) - mesh.x(pnb);
			k = mesh.y(idx) - mesh.y(pnb);
			weight[i] = exp(-3.0 * (h*h+k*k)/sl2 );
			A[mi+0] = 0.0;
			A[mi+1] = 1.0 * weight[i];
			A[mi+2] = h   * weight[i];
			A[mi+3] = k   * weight[i];
			A[mi+4] = h*h * weight[i];
			A[mi+5] = k*k * weight[i];
			A[mi+6] = h*k * weight[i];
		}
		else
		{
			diagEntry = i;
			weight[i] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			A[mi+2] = 0.0;
			A[mi+3] = 0.0;
			A[mi+4] = 0.0;
			A[mi+5] = 0.0;
			A[mi+6] = 0.0;
		}
	}

	int n = 7;
	int m = NP;
	double work[m+n];
	int info, workl = m+n;
	int one = 1;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		entries(i) = b[i]*weight[i];
	}

	//~ // Check for diagonally dominance
	//~ double diagEnt = abs(entries(diagEntry));
	//~ double sumte = 0;
	//~ cout << "Entry: " << pnb << endl;
	//~ for (int i=0; i<NP; i++)
	//~ {
	//~ if ( i != diagEntry )
	//~ {
	//~ sumte += abs(entries(i));
	//~ cout << entries(i) << " ";
	//~ }
	//~ else
	//~ cout << " * " << entries(i) << " * ";
	//~ }
	//~ cout << endl;
	//~ cout << "LAPDiff: " << diagEnt << " > " << sumte << " = " << (sumte < diagEnt) << endl;



	// Test operator:
	double sum = 0;
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		sum += entries(i)*(mesh.x(idx)*mesh.x(idx)+mesh.y(idx)*mesh.y(idx));
	}
	if ( abs(sum - 4.0) > 1e-7 )
	{
		return -2;
	}

	return 0;


}



int partialLaplaceMatEntry(const pointSet& mesh, int pnb, Array& active, Array& entries)
{
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[7*NP];
	double weight[NP];
	int arrayEntry[NP];

	if (NP < 7)
		throw mfpmExcept(20);

	b[0] = -10.0 / pow(mesh.getSL(pnb),2);
	b[1] = 0.0;
	b[2] = 0.0;
	b[3] = 0.0;
	b[4] = 2.0;
	b[5] = 2.0;
	b[6] = 0.0;

	int diagEntry;
	int idx, mi, cnt;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	cnt = 0;
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = cnt*7;
		if (idx!=pnb)
		{
			if ( active(i) > 0.5 )
			{
				h = mesh.x(idx) - mesh.x(pnb);
				k = mesh.y(idx) - mesh.y(pnb);
				weight[cnt] = exp(-3.0 * (h*h+k*k)/sl2 );
				A[mi+0] = 0.0;
				A[mi+1] = 1.0 * weight[cnt];
				A[mi+2] = h   * weight[cnt];
				A[mi+3] = k   * weight[cnt];
				A[mi+4] = h*h * weight[cnt];
				A[mi+5] = k*k * weight[cnt];
				A[mi+6] = h*k * weight[cnt];
				arrayEntry[cnt] = i;
				cnt++;
			}
		}
		else
		{
			diagEntry = i;
			weight[cnt] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			A[mi+2] = 0.0;
			A[mi+3] = 0.0;
			A[mi+4] = 0.0;
			A[mi+5] = 0.0;
			A[mi+6] = 0.0;
			arrayEntry[cnt] = i;
			cnt++;
		}
	}

	if (cnt < 7)
		throw mfpmExcept(20);


	int n = 7;
	int m = cnt;
	double work[2*(m+n)];
	int info, workl = 2*(m+n);
	int one = 1;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<cnt; i++)
	{
		entries(arrayEntry[i]) = b[i]*weight[i];
	}


	// Test operator:
	//~ double sum = 0;
	//~ for (int i=0; i<NP; i++)
	//~ {
	//~ idx = mesh.getIDX(pnb,i);
	//~ sum += entries(i)*(mesh.x(idx)*mesh.x(idx)+mesh.y(idx)*mesh.y(idx));
	//~ }
	//~ if ( abs(sum - 4.0) > 1e-7 )
	//~ {
	//~ return -2;
	//~ }

	return 0;


}



int laplaceOrder4MatEntry(const pointSet& mesh, int pnb, Array& entries)
{
	int nbEq = 15;
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[nbEq*NP];
	double weight[NP];

	if (NP < nbEq)
		throw mfpmExcept(20);

	b[0] = 0.0;		// const
	b[1] = 0.0;		// h
	b[2] = 0.0;		// k
	b[3] = 2.0;		// hh
	b[4] = 2.0;		// kk
	b[5] = 0.0;		// hk

	b[6] = 0.0;		// hhh
	b[7] = 0.0;		// hhk
	b[8] = 0.0;		// hkk
	b[9] = 0.0;		// kkk

	b[10] = 0.0;	// hhhh
	b[11] = 0.0;	// kkkk
	b[12] = 0.0;	// hhhk
	b[13] = 0.0;	// hhkk
	b[14] = 0.0;	// hkkk


	double maxLength = 0;
	for (int i=0; i<NP; i++)
	{
		int idx = mesh.getIDX(pnb,i);
		double dist = mesh.P(pnb).dist(mesh.P(idx));
		if (dist>maxLength)
			maxLength=dist;
	}
	b[3] /= maxLength*maxLength;
	b[4] /= maxLength*maxLength;


	int diagEntry;
	int idx, mi;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = nbEq*i;

		h = mesh.x(idx) - mesh.x(pnb);
		k = mesh.y(idx) - mesh.y(pnb);
		weight[i] = exp(-2.0 * (h*h+k*k)/sl2 );

		h /= maxLength;
		k /= maxLength;

		A[mi+0] = 1.0 * weight[i];
		A[mi+1] = h   * weight[i];
		A[mi+2] = k   * weight[i];
		A[mi+3] = h*h * weight[i];
		A[mi+4] = k*k * weight[i];
		A[mi+5] = h*k * weight[i];

		double fac=1.0;
		A[mi+6] = fac*h*h*h * weight[i];
		A[mi+7] = fac*h*h*k * weight[i];
		A[mi+8] = fac*h*k*k * weight[i];
		A[mi+9] = fac*k*k*k * weight[i];

		A[mi+10] = h*h*h*h * weight[i];
		A[mi+11] = k*k*k*k * weight[i];
		A[mi+12] = h*h*h*k * weight[i];
		A[mi+13] = h*h*k*k * weight[i];
		A[mi+14] = h*k*k*k * weight[i];

	}



	int n = nbEq;
	int lda = nbEq;
	int m = NP;
	int jpv[m+n];
	int lwork = 3*(m+n);
	double work[lwork];
	int info;
	int one = 1;
	int rank;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &lda, b, &m, work, &lwork, &info);


	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		entries(i) = b[i]*weight[i];
	}

	// Test operator:
	double sum = 0;
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		sum += entries(i)*(mesh.x(idx)*mesh.x(idx)+mesh.y(idx)*mesh.y(idx));
	}
	if ( abs(sum - 4.0) > 1e-7 )
	{
		//cout << pnb << "   " << sum << endl;
		return -2;
	}


	return 0;


}



int biLaplaceMatEntry(const pointSet& mesh, int pnb, Array& entries)
{
	int nbEq = 15;
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[nbEq*NP];
	double weight[NP];

	if (NP < nbEq)
		throw mfpmExcept(20);

	b[0] = 0.0;		// const
	b[1] = 0.0;		// h
	b[2] = 0.0;		// k
	b[3] = 0.0;		// hh
	b[4] = 0.0;		// kk
	b[5] = 0.0;		// hk

	b[6] = 0.0;		// hhh
	b[7] = 0.0;		// hhk
	b[8] = 0.0;		// hkk
	b[9] = 0.0;		// kkk

	b[10] = 24.0;	// hhhh
	b[11] = 24.0;	// kkkk
	b[12] = 0.0;	// hhhk
	b[13] = 8.0;	// hhkk
	b[14] = 0.0;	// hkkk

	int diagEntry;
	int idx, mi;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = nbEq*i;

		h = mesh.x(idx) - mesh.x(pnb);
		k = mesh.y(idx) - mesh.y(pnb);
		weight[i] = exp(-2.0 * (h*h+k*k)/sl2 );

		A[mi+0] = 1.0 * weight[i];
		A[mi+1] = h   * weight[i];
		A[mi+2] = k   * weight[i];
		A[mi+3] = h*h * weight[i];
		A[mi+4] = k*k * weight[i];
		A[mi+5] = h*k * weight[i];

		double fac=1.0;
		A[mi+6] = fac*h*h*h * weight[i];
		A[mi+7] = fac*h*h*k * weight[i];
		A[mi+8] = fac*h*k*k * weight[i];
		A[mi+9] = fac*k*k*k * weight[i];

		A[mi+10] = h*h*h*h * weight[i];
		A[mi+11] = k*k*k*k * weight[i];
		A[mi+12] = h*h*h*k * weight[i];
		A[mi+13] = h*h*k*k * weight[i];
		A[mi+14] = h*k*k*k * weight[i];

	}



	int n = nbEq;
	int lda = nbEq;
	int m = NP;
	int jpv[m+n];
	int lwork = 3*(m+n);
	double work[lwork];
	int info;
	int one = 1;
	int rank;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &lda, b, &m, work, &lwork, &info);


	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		entries(i) = b[i]*weight[i];
	}

	// Test operator:
	double sum = 0;
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		sum += entries(i)*(pow(mesh.x(idx),4)+pow(mesh.y(idx),4));
	}
	if ( abs(sum - 48.0) > 1e-7 )
	{
		cout << pnb << "   " << sum << endl;
		return -2;
	}


	return 0;


}



int gradMatEntry(const pointSet& mesh, int pnb, Array& gradX, Array& gradY, int order)
{
	int orderFac;
	int NP = mesh.getNIDX(pnb);

	if (order==1 & NP>=4)
		orderFac = 4;
	else if (order==2 & NP>=7)
		orderFac = 7;
	else if (order==3 & NP>=11)
		orderFac = 11;
	else if (order==4 & NP>=16)
		orderFac = 16;
	else
		throw mfpmExcept(20);

	double b[2*NP];
	double A[16*NP];
	double weight[NP];

	if (NP < 4)
		throw mfpmExcept(20);

	double sl = mesh.getSL(pnb);
	double sl2 = pow(sl,2);
	double alpha = sl;

	// GradX
	b[0] = 0.0;
	b[1] = 0.0;
	b[2] = 1.0 / alpha;
	for (int i=3; i<orderFac; i++)
		b[i] = 0.0;

	// GradY
	b[NP+0] = 0.0;
	b[NP+1] = 0.0;
	b[NP+2] = 0.0;
	b[NP+3] = 1.0 / alpha;
	for (int i=4; i<orderFac; i++)
		b[NP+i] = 0.0;

	int diagEntry;
	int idx, mi;
	double h, k;
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = orderFac*i;
		if (idx!=pnb)
		{
			h = mesh.x(idx) - mesh.x(pnb);
			k = mesh.y(idx) - mesh.y(pnb);
			weight[i] = exp(-2.0 * (h*h+k*k)/sl2 );
			h /= alpha;
			k /= alpha;
			A[mi+0] = 0.0;
			A[mi+1] = 1.0 * weight[i];
			A[mi+2] = h   * weight[i];
			A[mi+3] = k   * weight[i];
			if (orderFac>4)
			{
				A[mi+4] = h*h * weight[i];
				A[mi+5] = k*k * weight[i];
				A[mi+6] = h*k * weight[i];
			}
			if (orderFac>7)
			{
				A[mi+7]  = h*h*h *weight[i];
				A[mi+8]  = h*h*k *weight[i];
				A[mi+9]  = h*k*k *weight[i];
				A[mi+10] = k*k*k *weight[i];
			}
			if (orderFac>11)
			{
				A[mi+11] = h*h*h*h *weight[i];
				A[mi+12] = h*h*h*k *weight[i];
				A[mi+13] = h*h*k*k *weight[i];
				A[mi+14] = h*k*k*k *weight[i];
				A[mi+15] = k*k*k*k *weight[i];
			}


		}
		else
		{
			diagEntry = i;
			weight[i] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			for (int i=2; i<orderFac; i++)
				A[mi+i] = 0.0;

		}
	}

	int n = orderFac;
	int m = NP;
	int workl = 3*(m+n);
	double work[workl];
	int info;
	int nrhs = 2;
	char TN[] = "N";
	dgels_(TN, &n, &m, &nrhs, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		gradX(i) = b[i]*weight[i];
		gradY(i) = b[NP+i]*weight[i];
	}

	// Test operator:
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		sum1 += gradX(i)*(mesh.x(idx)+mesh.y(idx));
		sum2 += gradY(i)*(mesh.x(idx)+mesh.y(idx));
	}
	if ( (sum1 - 1.0) > 1e-7 )
	{
		return -1;
	}
	if ( (sum2 - 1.0) > 1e-7 )
	{
		return -2;
	}


	return 0;



}



int neumannMatEntry(const pointSet& mesh, int pnb, Array& entries, bool forceLowOrder)
{
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[6*NP];
	double weight[NP];

	if (NP < 4)
		throw mfpmExcept(20);

	double sl2 = pow(mesh.getSL(pnb),2);

	b[0] = 3.0 / mesh.getSL(pnb);
	b[1] = 0.0;
	b[2] = mesh.getNormal(0,pnb);
	b[3] = mesh.getNormal(1,pnb);
	b[4] = 0.0;
	b[5] = 0.0;
	int orderFac = 6;
	if (NP<6 | forceLowOrder)
		orderFac = 4;

	int idx, mi;
	int diagEntry=0;
	double h, k;
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = orderFac*i;
		if (idx!=pnb)
		{
			h = mesh.x(idx) - mesh.x(pnb);
			k = mesh.y(idx) - mesh.y(pnb);
			weight[i] = exp(-2.0 * (h*h+k*k)/sl2 );
			A[mi+0] = 0.0;
			A[mi+1] = 1.0 * weight[i];
			A[mi+2] = h   * weight[i];
			A[mi+3] = k   * weight[i];
			if (orderFac>4)
			{
				A[mi+4] = h*h * weight[i];
				A[mi+5] = k*k * weight[i];
			}
		}
		else
		{
			diagEntry = i;
			weight[i] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			A[mi+2] = 0.0;
			A[mi+3] = 0.0;
			if (orderFac>4)
			{
				A[mi+4] = 0.0;
				A[mi+5] = 0.0;
			}
		}
	}

	int n = orderFac;
	int m = NP;
	double work[m+n];
	int info, workl = m+n;
	int one = 1;
	char TN[] ="N";
	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		entries(i) = b[i]*weight[i];
	}


	// Test operator:
	double tmp = entries(diagEntry);
	for (int i=0; i<NP; i++)
	{
		if (abs(entries(i))>tmp)
		{
			//cout << "Invalid Neumann matrix entries!!" << endl;
			//return -1;
		}
	}
	double sum = 0;
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		sum += entries(i)*(mesh.x(idx)+mesh.y(idx));
	}
	if ( (sum - mesh.getNormal(0,pnb) - mesh.getNormal(1,pnb)) > 1e-7 )
	{
		//cout << "Invalid Neumann matrix entries!!" << endl;
		return -2;
	}

	return 0;

}



void approximatedNormals(const pointSet& mesh, vector<int> neighb, int pnb, double& nx_out, double& ny_out)
{

	if (mesh.getBND(pnb)==0)
		throw mfpmExcept(0);

	int NP = neighb.size();
	double sl2 = mesh.getSL(pnb);


	// New approach
	double sumX = 0, sumY = 0;
	double nx, ny;
	int count = 0;
	for (int i=0; i<NP; i++)
	{
		if (pnb!=neighb[i] )
		{
			nx = (mesh.y(pnb)-mesh.y(neighb[i]));
			ny = -(mesh.x(pnb)-mesh.x(neighb[i]));
			double weight = exp(-5.0* (nx*nx+ny*ny) / sl2);

			double length = sqrt(nx*nx+ny*ny);
			nx /= length;
			ny /= length;
			if ( nx*mesh.getNormal(0,pnb) + ny*mesh.getNormal(1,pnb) < 0)
			{
				nx *= -1.0;
				ny *= -1.0;
			}
			sumX += nx*weight;
			sumY += ny*weight;
		}
	}
	sumX /= NP-1;
	sumY /= NP-1;
	double length = sqrt(sumX*sumX+sumY*sumY);
	if (length==0)
		throw mfpmExcept(1000);
	nx_out = sumX / length;
	ny_out = sumY / length;


}



int interpolationConstants(const scalarField& f, int pos, double* con, bool useInner)
{
	const pointSet* mesh = f.getMesh();
	const point& P = mesh->P(pos);
	vector<int> neighb;
	//mesh->findPossibleNeighbours(P, neighb, useInner);
	for (int i=0; i<mesh->getNIDX(pos); i++)
		neighb.push_back(mesh->getIDX(pos,i));
	int m = neighb.size();
	double rhs[m];
	int mi;
	double A[15*m];
	double weight;
	double h2 = pow((mesh->getMaxSL()+mesh->getMinSL())/2,2);
	int orderFac = 15;
	if (m>14)
		orderFac = 15;
	else if (m>9)
		orderFac = 10;
	else if (m>5)
		orderFac = 6;
	else if (m>2)
		orderFac = 3;
	else
		throw mfpmExcept(20);



	for (int i=0; i<m; i++)
	{
		const point& Ptmp = mesh->P(neighb[i]);
		mi = orderFac*i;
		weight = exp(-2.00*P.dist2(Ptmp) / h2);
		A[mi+0] = weight*1.0;
		A[mi+1] = weight*Ptmp.X;
		A[mi+2] = weight*Ptmp.Y;
		if (orderFac>3)
		{
			A[mi+3] = weight*Ptmp.X*Ptmp.X;
			A[mi+4] = weight*Ptmp.X*Ptmp.Y;
			A[mi+5] = weight*Ptmp.Y*Ptmp.Y;
		}
		if (orderFac>6)
		{
			A[mi+6] = weight*Ptmp.X*Ptmp.X*Ptmp.X;
			A[mi+7] = weight*Ptmp.X*Ptmp.X*Ptmp.Y;
			A[mi+8] = weight*Ptmp.X*Ptmp.Y*Ptmp.Y;
			A[mi+9] = weight*Ptmp.Y*Ptmp.Y*Ptmp.Y;
		}
		if (orderFac>10)
		{
			A[mi+10] = weight*Ptmp.X*Ptmp.X*Ptmp.X*Ptmp.X;
			A[mi+11] = weight*Ptmp.X*Ptmp.X*Ptmp.X*Ptmp.Y;
			A[mi+12] = weight*Ptmp.X*Ptmp.X*Ptmp.Y*Ptmp.Y;
			A[mi+13] = weight*Ptmp.X*Ptmp.Y*Ptmp.Y*Ptmp.Y;
			A[mi+14] = weight*Ptmp.Y*Ptmp.Y*Ptmp.Y*Ptmp.Y;
		}
		rhs[i]  = weight*f(neighb[i]);
	}

	int n = orderFac;
	double work[m+n];
	int info, workl = m+n;
	int one = 1;
	char TN[] = "T";
	dgels_(TN, &n, &m, &one, A, &n, rhs, &m, work, &workl, &info);

	for (int i=0; i<orderFac; i++)
		con[i] = rhs[i];

	if (orderFac == 3)
		return 1;
	else if (orderFac == 6)
		return 2;
	else if (orderFac == 10)
		return 3;
	else if (orderFac > 10)
		return 4;

	return -1;


}





int partialNeumannMatEntry(const pointSet& mesh, int pnb, const Array& active, Array& entries, bool forceLowOrder)
{
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[6*NP];
	double weight[NP];
	int arrayEntry[NP];

	if (NP < 4)
		throw mfpmExcept(20);

	double sl2 = pow(mesh.getSL(pnb),2);

	b[0] = 0*3.0 / mesh.getSL(pnb);
	b[1] = 0.0;
	b[2] = mesh.getNormal(0,pnb);
	b[3] = mesh.getNormal(1,pnb);
	b[4] = 0.0;
	b[5] = 0.0;
	b[6] = 0.0;
	int orderFac = 6;
	if (NP<6 | forceLowOrder)
		orderFac = 4;

	int idx, mi;
	int diagEntry=0;
	double h, k;
	int cnt = 0;
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = orderFac*cnt;
		if (idx!=pnb)
		{
			if (active(i) > 0.5)
			{
				h = mesh.x(idx) - mesh.x(pnb);
				k = mesh.y(idx) - mesh.y(pnb);
				weight[cnt] = exp(-2.0 * (h*h+k*k)/sl2 );
				A[mi+0] = 0.0;
				A[mi+1] = 1.0 * weight[cnt];
				A[mi+2] = h   * weight[cnt];
				A[mi+3] = k   * weight[cnt];
				if (orderFac>4)
				{
					A[mi+4] = h*h * weight[cnt];
					A[mi+5] = h*k * weight[cnt];
					A[mi+6] = k*k * weight[cnt];
				}
				arrayEntry[cnt] = i;
				cnt++;
			}
		}
		else
		{
			weight[cnt] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			A[mi+2] = 0.0;
			A[mi+3] = 0.0;
			if (orderFac>4)
			{
				A[mi+4] = 0.0;
				A[mi+5] = 0.0;
			}
			arrayEntry[cnt] = i;
			cnt++;
		}
	}

	int n = orderFac;
	int m = cnt;
	double work[2*(m+n)];
	int info, workl = 2*(m+n);
	int one = 1;
	char TN[] ="N";


	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	entries = 0.0;
	for (int i=0; i<cnt; i++)
	{
		entries(arrayEntry[i]) = b[i]*weight[i];
	}


	//~ // Test operator:
	//~ double tmp = entries(diagEntry);
	//~ for (int i=0; i<NP; i++)
	//~ {
	//~ if (abs(entries(i))>tmp)
	//~ {
	//~ //cout << "Invalid Neumann matrix entries!!" << endl;
	//~ //return -1;
	//~ }
	//~ }
	//~ double sum = 0;
	//~ for (int i=0; i<NP; i++)
	//~ {
	//~ idx = mesh.getIDX(pnb,i);
	//~ sum += entries(i)*(mesh.x(idx)+mesh.y(idx));
	//~ }
	//~ if ( (sum - mesh.getNormal(0,pnb) - mesh.getNormal(1,pnb)) > 1e-7 )
	//~ {
	//~ //cout << "Invalid Neumann matrix entries!!" << endl;
	//~ return -2;
	//~ }

	return 0;

}




int partialGradMatEntry(const pointSet& mesh, int pnb, const Array& active, Array& gradX, Array& gradY, int order)
{
	int orderFac;
	int NP = mesh.getNIDX(pnb);
	int cntNP = 0;
	for (int i=0; i<NP; i++)
	{
		if (active(i)>0.5)
			cntNP++;
	}

	if (order==1 & cntNP>=4)
		orderFac = 4;
	else if (order==2 & cntNP>=7)
		orderFac = 7;
	else if (order==3 & cntNP>=11)
		orderFac = 11;
	else if (order==4 & cntNP>=16)
		orderFac = 16;
	else
		throw mfpmExcept(20);

	double b[2*NP];
	double A[16*NP];
	double weight[NP];
	int entryArray[NP];

	if (NP < 4)
		throw mfpmExcept(20);

	double sl = mesh.getSL(pnb);
	double sl2 = pow(sl,2);
	double alpha = sl;


	int diagEntry;
	int idx, mi, cnt;
	double h, k;
	cnt = 0;
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = orderFac*cnt;
		if (idx!=pnb)
		{
			if (active(i) > 0.5)
			{
				h = mesh.x(idx) - mesh.x(pnb);
				k = mesh.y(idx) - mesh.y(pnb);
				weight[cnt] = exp(-2.0 * (h*h+k*k)/sl2 );
				h /= alpha;
				k /= alpha;
				A[mi+0] = 0.0;
				A[mi+1] = 1.0 * weight[cnt];
				A[mi+2] = h   * weight[cnt];
				A[mi+3] = k   * weight[cnt];
				if (orderFac>4)
				{
					A[mi+4] = h*h * weight[cnt];
					A[mi+5] = k*k * weight[cnt];
					A[mi+6] = h*k * weight[cnt];
				}
				entryArray[cnt] = i;
				cnt++;
			}
		}
		else
		{
			diagEntry = i;
			weight[i] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			for (int j=2; j<orderFac; j++)
				A[mi+j] = 0.0;
			entryArray[cnt] = i;
			cnt++;
		}
	}

	// GradX
	b[0] = 0.0;
	b[1] = 0.0;
	b[2] = 1.0 / alpha;
	for (int i=3; i<orderFac; i++)
		b[i] = 0.0;

	// GradY
	b[cnt+0] = 0.0;
	b[cnt+1] = 0.0;
	b[cnt+2] = 0.0;
	b[cnt+3] = 1.0 / alpha;
	for (int i=4; i<orderFac; i++)
		b[cnt+i] = 0.0;

	int n = orderFac;
	int m = cnt;
	int workl = 3*(m+n);
	double work[workl];
	int info;
	int nrhs = 2;
	char TN[] = "N";
	dgels_(TN, &n, &m, &nrhs, A, &n, b, &m, work, &workl, &info);
	cout << "NP: " << cnt << endl;

	if (info!=0)
	{
		cout <<  "Error solving least squares: partialGradMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	gradX = 0;
	gradY = 0;
	for (int i=0; i<cnt; i++)
	{
		gradX(entryArray[i]) = b[i]*weight[i];
		gradY(entryArray[i]) = b[cnt+i]*weight[i];
	}

	return 0;



}






int partialApproximationMatEntry(const pointSet& mesh, int pnb, Array& active, Array& entries)
{
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[7*NP];
	double weight[NP];
	int arrayEntry[NP];

	if (NP < 5)
		throw mfpmExcept(20);

	int orderFac = 6;

	b[0] = 1.0;
	b[1] = 0.0;
	b[2] = 0.0;
	if (orderFac > 3)
	{
		b[3] = 0.0;
		b[4] = 0.0;
		b[5] = 0.0;
	}

	int diagEntry;
	int idx, mi, cnt;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	cnt = 0;
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = cnt*orderFac;
		if ( active(i) > 0.5 )
		{
			h = mesh.x(idx) - mesh.x(pnb);
			k = mesh.y(idx) - mesh.y(pnb);
			weight[cnt] = exp(-3.0 * (h*h+k*k)/sl2 );
			A[mi+0] = 1.0 * weight[cnt];
			A[mi+1] = h   * weight[cnt];
			A[mi+2] = k   * weight[cnt];
			if (orderFac>3)
			{
				A[mi+3] = h*h * weight[cnt];
				A[mi+4] = h*k   * weight[cnt];
				A[mi+5] = k*k   * weight[cnt];
			}
			arrayEntry[cnt] = i;
			cnt++;
		}
	}

	if (cnt < 4)
		throw mfpmExcept(20);


	int n = orderFac;
	int m = cnt;
	double work[2*(m+n)];
	int info, workl = 2*(m+n);
	int one = 1;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<cnt; i++)
	{
		entries(arrayEntry[i]) = b[i]*weight[i];
	}


	return 0;


}





int divLambdaGradMatEntry(const pointSet& mesh, int pnb, double lambda, double DxLambda, double DyLambda, const Array& larr, Array& entries)
{
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[7*NP];
	double weight[NP];

	if (NP < 7)
		throw mfpmExcept(20);

	double nx = DxLambda;
	double ny = DyLambda;
	double nabs = sqrt(nx*nx+ny*ny);
	if ( nabs > 1e-5 )
	{
		nx /= nabs;
		ny /= nabs;
	}
	else
	{
		nx = 1.0;
		ny = 0.0;
	}
	double tx = -ny;
	double ty = nx;

	b[0] = -10 / pow(mesh.getSL(pnb),2);
	b[1] = 0.0;
	b[2] = DxLambda;
	b[3] = DyLambda;
	b[4] = 2.0*lambda;
	b[5] = 2.0*lambda;
	b[6] = 0.0;

	double alphaL = 3.0;
	double alphaR = 10.0;



	int diagEntry;
	int idx, mi;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = 7*i;
		if (idx!=pnb)
		{
			h = mesh.x(idx) - mesh.x(pnb);
			k = mesh.y(idx) - mesh.y(pnb);
			double hN = h*nx + k*ny;
			double hT = h*tx + k*ty;
			weight[i] = exp(-alphaL*hT*hT/sl2 );
			double alpha = alphaL;
			// TODO: Only for particular test case!!!!!
			if ( lambda < 0.1 )
				if ( hN > 0)
					alpha = alphaR;
				else
					alpha = alphaL;
			else if (hN > 0)
				alpha = alphaL;
			else
				alpha = alphaR;
			weight[i] *= exp(-alpha*hN*hN/sl2 );
			A[mi+0] = 0.0;
			A[mi+1] = 1.0 * weight[i];
			A[mi+2] = h   * weight[i];
			A[mi+3] = k   * weight[i];
			A[mi+4] = h*h * weight[i];
			A[mi+5] = k*k * weight[i];
			A[mi+6] = h*k * weight[i];
		}
		else
		{
			diagEntry = i;
			weight[i] = 1.0;
			A[mi+0] = 1.0;
			A[mi+1] = 1.0;
			A[mi+2] = 0.0;
			A[mi+3] = 0.0;
			A[mi+4] = 0.0;
			A[mi+5] = 0.0;
			A[mi+6] = 0.0;
		}
	}

	int n = 7;
	int m = NP;
	double work[m+n];
	int info, workl = m+n;
	int one = 1;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		entries(i) = b[i]*weight[i];
	}

	// Check for diagonally dominance
	double diagEnt = abs(entries(diagEntry));
	double sum = 0;
	cout << "Entry: " << pnb << endl;
	for (int i=0; i<NP; i++)
	{
		if ( i != diagEntry )
		{
			sum += abs(entries(i));
			cout << entries(i) << " ";
		}
		else
			cout << " * " << entries(i) << " * ";
	}
	cout << endl;
	cout << "Diff: " << diagEnt << " > " << sum << " = " << (sum < diagEnt) << endl;


	return 0;


}




int approximationXYMatEntry(const pointSet& mesh, int pnb, double xp, double yp, Array& entries, int order)
{
	int orderFac = 3;
	int NP = mesh.getNIDX(pnb);
	double b[NP];
	double A[7*NP];
	double weight[NP];
	int arrayEntry[NP];

	if (order == 1)
		orderFac = 3;
	else
		orderFac = 6;

	if (NP < orderFac)
		throw mfpmExcept(20);

	b[0] = 1.0;
	b[1] = 0.0;
	b[2] = 0.0;

	b[3] = 0.0;
	b[4] = 0.0;
	b[5] = 0.0;


	int diagEntry;
	int idx, mi;
	double h, k;
	double sl2 = pow(mesh.getSL(pnb),2);
	// Store matrix in Fortran format
	for (int i=0; i<NP; i++)
	{
		idx = mesh.getIDX(pnb,i);
		mi = i*orderFac;
		h = mesh.x(idx) - xp;
		k = mesh.y(idx) - yp;
		weight[i] = exp(-3.0 * (h*h+k*k)/sl2 );
		A[mi+0] = 1.0 * weight[i];
		A[mi+1] = h   * weight[i];
		A[mi+2] = k   * weight[i];
		if (orderFac > 3)
		{
			A[mi+3] = h*h   * weight[i];
			A[mi+4] = k*k   * weight[i];
			A[mi+5] = h*k   * weight[i];
		}

	}



	int n = orderFac;
	int m = NP;
	double work[2*(m+n)];
	int info, workl = 2*(m+n);
	int one = 1;
	char TN[] = "N";
	dgels_(TN, &n, &m, &one, A, &n, b, &m, work, &workl, &info);

	if (info!=0)
	{
		cout <<  "Error solving least squares: laplaceMatEntry!" << endl;
		throw mfpmExcept(27);
	}

	for (int i=0; i<NP; i++)
	{
		entries(i) = b[i]*weight[i];
	}


	return 0;


}







