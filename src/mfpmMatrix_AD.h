
#ifndef MFPMMATRIX_AD_H
#define MFPMMATRIX_AD_H

#include "pointSet.h"
#include "Array_AD.h"
#include "scalarField.h"
#include "primitive_operators.h"

#include <adolc/adolc.h>
#include "mfpmMatrix.h"

class mfpmMatrix_AD
{
private:
	const pointSet* mesh;
	adouble** data;
	Array_AD* rhs;
	int n;

public:
	mfpmMatrix_AD(const pointSet& mesh);
	mfpmMatrix_AD(const mfpmMatrix_AD& mat);
	mfpmMatrix_AD(const mfpmMatrix& mat);
	mfpmMatrix_AD();
	~mfpmMatrix_AD();

	void clearBnd();
	void setDirichlet(int bnd, double val);
	void setDirichlet(int bnd, const Array_AD& val);
	void setNeumann(int bnd, double val);
	void setNeumann(int bnd, const Array_AD& val);
	void setRobin(int bnd, double eps, double val);
	void setRobin(int bnd, double eps, const Array_AD& val);
	void addRhs(int i, adouble& val);
	void setRhs(int i, adouble& val);


	adouble operator()(int i, int j) const;
	void operator=(const mfpmMatrix_AD& mat);
	void operator=(const mfpmMatrix& mat);
	// Multiplying
	Array_AD operator*(const Array_AD& x) const;
	Array_AD operator*(const scalarField& x) const;
	mfpmMatrix_AD operator*(const adouble& d) const;
	// Add
	mfpmMatrix_AD operator+(const Array_AD& V) const;
	mfpmMatrix_AD operator+(const adouble& d) const;
	mfpmMatrix_AD operator+(const mfpmMatrix_AD& mat) const;
	// Sub
	mfpmMatrix_AD operator-(const mfpmMatrix_AD& mat) const;
	mfpmMatrix_AD operator-(const adouble& d) const;
	mfpmMatrix_AD operator-(const Array_AD& V) const;

	// Inline
	adouble getData(int i, int j) const {
		return data[i][j];
	};
	adouble** getData() const {
		return data;
	};
	Array_AD* getRhs() const {
		return rhs;
	};
	adouble getRhs(int i) const {
		return rhs->value(i);
	};
	const pointSet* getMesh() const {
		return mesh;
	};
	int getSize() const {
		return n;
	};

	// Friends
	friend mfpmMatrix_AD operator*(const adouble& d, const mfpmMatrix_AD& mat);
	friend mfpmMatrix_AD operator*(const Array_AD& d, const mfpmMatrix_AD& mat);
	friend mfpmMatrix_AD operator+(const adouble& d, const mfpmMatrix_AD& mat);
	friend mfpmMatrix_AD operator+(const Array_AD& d, const mfpmMatrix_AD& mat);
	friend mfpmMatrix_AD operator-(const adouble& d, const mfpmMatrix_AD& mat);
	friend mfpmMatrix_AD operator-(const Array_AD& d, const mfpmMatrix_AD& mat);
	friend mfpmMatrix_AD operator==(const mfpmMatrix_AD& mat, const Array_AD& V);
	friend mfpmMatrix_AD operator==(const mfpmMatrix_AD& mat, const adouble& d);


	void setRow(int row, Array_AD val);
	void setDiag(int row, adouble& val);
	void addRow(int row, Array_AD val);
	Array_AD getRow(int row) const;
	void regularise(double factor = 1e-4);
	void normalise();




};



#endif
