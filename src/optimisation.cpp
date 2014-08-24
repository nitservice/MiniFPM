#include "optimisation.h"

#include <iostream>
#include <cmath>

#include "Array.h"

using namespace std;

BFGS::BFGS(int size, double initVal)
{
	n_ = size;
	A_ = new double[size*size];
	setIdentity(initVal);
}


BFGS::~BFGS()
{
	delete A_;
}


void BFGS::setIdentity(double val)
{
	for (int i=0; i<n_; i++)
		for (int j=0; j<n_; j++)
			setHess(i,j,0.0);
	for (int i=0; i<n_; i++)
		setHess(i,i,val);
}


void BFGS::setHess(int i, int j, double val)
{
	A_[i + n_*j] = val;
}

double BFGS::getHess(int i, int j) const
{
	return A_[i + n_*j];
}


void BFGS::rankOneUpdate(double* xold, double* xnew, double* gradOld, double* gradNew)
{
	double s[n_];
	double y[n_];
	double sTs = 0.0;
	// Copy values
	for (int i=0; i<n_; i++)
	{
		s[i] = xnew[i] - xold[i];
		y[i] = gradNew[i] - gradOld[i];
	}
	// Calc s^Ts
	for (int i=0; i<n_; i++)
		sTs += s[i]*s[i];
	// Perform update
	for (int i=0; i<n_; i++)
	{
		double Hs = 0.0;
		for (int k=0; k<n_; k++)
		{
			Hs += getHess(i,k)*s[k];
		}
		for (int j=0; j<n_; j++)
		{
			double Hij = getHess(i,j) + (y[i]-Hs)*s[j]/sTs;
			setHess(i,j,Hij);
		}
	}
}


void BFGS::rankOneUpdate(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew)
{
	double a1[n_];
	double a2[n_];
	double a3[n_];
	double a4[n_];

	for (int i=0; i<n_; i++)
	{
		a1[i] = xold(i);
		a2[i] = xnew(i);
		a3[i] = gradOld(i);
		a4[i] = gradNew(i);
	}
	rankOneUpdate(a1, a2, a3 ,a4);
}


void BFGS::rankTwoUpdate(double* xold, double* xnew, double* gradOld, double* gradNew)
{
	cout << "Update Hessian with BFGS...   " << flush;
	double s[n_];
	double y[n_];
	double sTs = 0.0;
	// Copy values
	for (int i=0; i<n_; i++)
	{
		s[i] = xnew[i] - xold[i];
		y[i] = gradNew[i] - gradOld[i];
	}

	double yTs = 0.0;
	for (int i=0; i<n_; i++)
		yTs += y[i]*s[i];

	double Hs[n_];
	for (int i=0; i<n_; i++)
	{
		Hs[i] = 0.0;
		for (int k=0; k<n_; k++)
		{
			Hs[i] += getHess(i,k)*s[k];
		}
	}

	if (yTs <= 0)
	{
		cout << "Hessian reset!" << endl;
		setIdentity();
		return;
	}

	double sTHs = 0.0;
	for (int i=0; i<n_; i++)
		sTHs += s[i]*Hs[i];

	for (int i=0; i<n_; i++)
	{
		for (int j=0; j<n_; j++)
		{
			double hij = getHess(i,j) + y[i]*y[j]/yTs - (Hs[i]*Hs[j])/sTHs;
			setHess(i, j, hij);
		}
	}

	cout << "done!" << endl;
}



void BFGS::rankTwoUpdateMBFGS(double* xold, double* xnew, double* gradOld, double* gradNew)
{
	cout << "Update Hessian with MBFGS...   " << flush;
	double s[n_];
	double y[n_];
	double sTs = 0.0;

	// norms and step-size for mbfgs
	double gradNorm = 0;
	double sNorm2   = 0;
	for (int i=0; i<n_; i++)
	{
		gradNorm += gradNew[i]*gradNew[i];
		sNorm2   += pow( xnew[i]-xold[i] , 2);
	}
	gradNorm = sqrt(gradNorm);

	double tk = 0;
	for (int i=0; i<n_; i++)
	{
		tk -= (gradNew[i] - gradOld[i]) * (xnew[i] - xold[i]);
	}
	tk = tk / sNorm2;
	tk = max(tk,0.0);
	tk = 1 + tk;


	// Copy values
	for (int i=0; i<n_; i++)
	{
		s[i] = xnew[i] - xold[i];
		y[i] = gradNew[i] - gradOld[i] + tk*gradNorm*s[i];
	}

	double yTs = 0.0;
	for (int i=0; i<n_; i++)
		yTs += y[i]*s[i];

	double Hs[n_];
	for (int i=0; i<n_; i++)
	{
		Hs[i] = 0.0;
		for (int k=0; k<n_; k++)
		{
			Hs[i] += getHess(i,k)*s[k];
		}
	}

	if (yTs <= 0)
	{
		cout << "Hessian reset!" << endl;
		setIdentity();
		return;
	}

	double sTHs = 0.0;
	for (int i=0; i<n_; i++)
		sTHs += s[i]*Hs[i];

	for (int i=0; i<n_; i++)
	{
		for (int j=0; j<n_; j++)
		{
			double hij = getHess(i,j) + y[i]*y[j]/yTs - (Hs[i]*Hs[j])/sTHs;
			setHess(i, j, hij);
		}
	}

	cout << "done!" << endl;
}




void BFGS::rankTwoUpdateDamped(double* xold, double* xnew, double* gradOld, double* gradNew)
{
	cout << "Update Hessian with damped BFGS...   " << flush;
	double s[n_];
	double y[n_];
	double sTs = 0.0;

	// norms and step-size for mbfgs
	double gradNorm = 0;
	double sNorm2   = 0;
	for (int i=0; i<n_; i++)
	{
		gradNorm += gradNew[i]*gradNew[i];
		sNorm2   += pow( xnew[i]-xold[i] , 2);
	}
	gradNorm = sqrt(gradNorm);


	// Copy values
	for (int i=0; i<n_; i++)
	{
		s[i] = xnew[i] - xold[i];
		y[i] = gradNew[i] - gradOld[i];
	}

	double Hs[n_];
	for (int i=0; i<n_; i++)
	{
		Hs[i] = 0.0;
		for (int k=0; k<n_; k++)
		{
			Hs[i] += getHess(i,k)*s[k];
		}
	}

	double yTs = 0.0;
	for (int i=0; i<n_; i++)
		yTs += y[i]*s[i];

	double sBs = 0.0;
	for (int i=0; i<n_; i++)
	{
		sBs += s[i]*Hs[i];
	}

	double theta;
	if ( yTs >= 0.2*sBs )
		theta = 1.0;
	else
		theta = 0.8*sBs / (sBs-yTs);

	double rk[n_];
	for (int i=0; i<n_; i++)
	{
		rk[i] = theta*y[i] + (1-theta)*Hs[i];
	}


	if (yTs <= 0)
	{
		cout << "Hessian reset!" << endl;
		setIdentity();
		return;
	}

	double sTr = 0.0;
	for (int i=0; i<n_; i++)
	{
		sTr += s[i]*rk[i];
	}


	for (int i=0; i<n_; i++)
	{
		for (int j=0; j<n_; j++)
		{
			double hij = getHess(i,j) + rk[i]*rk[j]/(sTr) - (Hs[i]*Hs[j])/sBs;
			setHess(i, j, hij);
		}
	}
	cout << "done!" << endl;
}







void BFGS::rankTwoUpdate(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew)
{
	double a1[n_];
	double a2[n_];
	double a3[n_];
	double a4[n_];

	for (int i=0; i<n_; i++)
	{
		a1[i] = xold(i);
		a2[i] = xnew(i);
		a3[i] = gradOld(i);
		a4[i] = gradNew(i);
	}
	rankTwoUpdate(a1, a2, a3 ,a4);
}


void BFGS::rankTwoUpdateMBFGS(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew)
{
	double a1[n_];
	double a2[n_];
	double a3[n_];
	double a4[n_];

	for (int i=0; i<n_; i++)
	{
		a1[i] = xold(i);
		a2[i] = xnew(i);
		a3[i] = gradOld(i);
		a4[i] = gradNew(i);
	}
	rankTwoUpdateMBFGS(a1, a2, a3 ,a4);
}

void BFGS::rankTwoUpdateDamped(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew)
{
	double a1[n_];
	double a2[n_];
	double a3[n_];
	double a4[n_];

	for (int i=0; i<n_; i++)
	{
		a1[i] = xold(i);
		a2[i] = xnew(i);
		a3[i] = gradOld(i);
		a4[i] = gradNew(i);
	}
	rankTwoUpdateDamped(a1, a2, a3 ,a4);
}







void BFGS::descent(double* grad, double* desc)
{
	int pivot[n_];
	int info;
	int nrhs = 1;
	int n = n_;
	double gradTmp[n_];
	double* Atmp = new double[n_*n_];
	for (int i=0; i<n_; i++)
		gradTmp[i] = grad[i];
	for (int i=0; i<n_*n_; i++)
		Atmp[i] = A_[i];
	cout << "Calculating descent direction...   " << flush;
	dgesv_(&n, &nrhs, Atmp, &n, pivot, gradTmp, &n, &info);
	delete Atmp;
	if (info != 0)
		cout << "!!!!!!!!   Info LAPACK: " << info << " !!!!!" << endl;
	for (int i=0; i<n_; i++)
		desc[i] = gradTmp[i];

	cout << "done!" << endl;
}


void BFGS::descent(const Array& grad, Array& desc)
{
	int pivot[n_];
	int info;
	int nrhs = 1;
	int n = n_;
	double gradTmp[n_];
	double* Atmp = new double[n_*n_];
	for (int i=0; i<n_; i++)
		gradTmp[i] = grad(i);
	for (int i=0; i<n_*n_; i++)
		Atmp[i] = A_[i];
	cout << "Calculating descent direction...   " << flush;
	dgesv_(&n, &nrhs, Atmp, &n, pivot, gradTmp, &n, &info);
	delete Atmp;
	if (info != 0)
		cout << "!!!!!!!!   Info LAPACK: " << info << " !!!!!" << endl;
	for (int i=0; i<n_; i++)
		desc(i) = gradTmp[i];

	cout << "done!" << endl;
}


