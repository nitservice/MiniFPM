#ifndef OPTIMISATION_H
#define OPTIMISATION_H

class Array;

extern "C" void dgesv_(int* n, int* nrhs, double *a, int* lda, int *ipivot, double *b, int* ldb, int *info);

class BFGS
{
private:
	double* A_;
	int n_;

public:
	BFGS(int size, double initVal = 1.0);
	~BFGS();

	void setIdentity(double val = 1.0);
	void setHess(int i, int j, double val);
	double getHess(int i, int j) const;

	void rankOneUpdate(double* xold, double* xnew, double* gradOld, double* gradNew);
	void rankOneUpdate(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew);

	void rankTwoUpdate(double* xold, double* xnew, double* gradOld, double* gradNew);
	void rankTwoUpdate(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew);

	void rankTwoUpdateMBFGS(double* xold, double* xnew, double* gradOld, double* gradNew);
	void rankTwoUpdateMBFGS(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew);

	void rankTwoUpdateDamped(double* xold, double* xnew, double* gradOld, double* gradNew);
	void rankTwoUpdateDamped(const Array& xold, const Array& xnew, const Array& gradOld, const Array& gradNew);

	void descent(double* grad, double* desc);
	void descent(const Array& grad, Array& desc);


};




#endif
