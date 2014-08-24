#include "matrixUmfpack.h"

#include "mfpmException.h"
#include "mfpmMatrix.h"

#include "pointSet.h"
#include "scalarField.h"

#include <umfpack.h>


matrixUmfpack::matrixUmfpack(int size, int nbNonZero)
{
	n  = size;
	nz = nbNonZero;
	col = new int[n+1];
	row = new int[nz];
	val = new double[nz];
	rhs = new double[n];
}

matrixUmfpack::matrixUmfpack()
{
	n  = -1;
	nz = -1;
	col = NULL;
	row = NULL;
	val = NULL;
	rhs = NULL;
}

matrixUmfpack::~matrixUmfpack()
{
	if (n>0)
	{
		delete[] col;
		delete[] row;
		delete[] val;
		delete[] rhs;
	}
}




void matrixUmfpack::setMatrix(const mfpmMatrix& mat)
{
	if (mat.getSize() != n)
		throw "Error matrixUmfpack:: setMatrix: sizes differ";

	const pointSet* mesh = mat.getMesh();
	double** data = mat.getData();
	int* Ti = new int[nz];
	int* Tj = new int[nz];
	double* Tx = new double[nz];

	int counter = 0;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = i;
			Tj[counter] = mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}
	umfpack_di_triplet_to_col(n, n, counter, Ti, Tj, Tx, col, row, val, (int*) NULL);

	for (int i=0; i<n; i++)
		rhs[i] = mat.getRhs(i);

	delete[] Ti;
	delete[] Tj;
	delete[] Tx;

}



void matrixUmfpack::setMatrix2x2(
    mfpmMatrix& M11, mfpmMatrix& M12,
    mfpmMatrix& M21, mfpmMatrix& M22 )
{
	if (   2*M11.getSize() != n
	        | 2*M12.getSize() != n
	        | 2*M21.getSize() != n
	        | 2*M22.getSize() != n )
	{
		throw "Error matrixUMFPACK:: setMatrix2x2: sizes differ";
	}

	int n = this->n / 2;

	int* Ti = new int[nz];
	int* Tj = new int[nz];
	double* Tx = new double[nz];
	const pointSet* mesh = M11.getMesh();
	double** data = M11.getData();

	int counter = 0;
	// Add Matrix 11
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = i;
			Tj[counter] = mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M12
	data = M12.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = i;
			Tj[counter] = n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M21
	data = M21.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = n+i;
			Tj[counter] = mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M22
	data = M22.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = n+i;
			Tj[counter] = n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	umfpack_di_triplet_to_col(this->n, this->n, counter, Ti, Tj, Tx, col, row, val, (int*) NULL);

	// Add RHS
	for (int i=0; i<n; i++)
	{
		rhs[i] = M11.getRhs(i) + M12.getRhs(i);
		rhs[n+i] = M21.getRhs(i) + M22.getRhs(i);
	}

	delete[] Ti;
	delete[] Tj;
	delete[] Tx;

}



void matrixUmfpack::setMatrix3x3(
    mfpmMatrix& M11, mfpmMatrix& M12, mfpmMatrix& M13,
    mfpmMatrix& M21, mfpmMatrix& M22, mfpmMatrix& M23,
    mfpmMatrix& M31, mfpmMatrix& M32, mfpmMatrix& M33 )
{
	if (   3*M11.getSize() != n
	        | 3*M12.getSize() != n
	        | 3*M13.getSize() != n
	        | 3*M21.getSize() != n
	        | 3*M22.getSize() != n
	        | 3*M23.getSize() != n
	        | 3*M31.getSize() != n
	        | 3*M32.getSize() != n
	        | 3*M33.getSize() != n)
	{
		throw "Error matrixUMFPACK:: setMatrix3x3: sizes differ";
	}

	int n = this->n / 3;

	int* Ti = new int[nz];
	int* Tj = new int[nz];
	double* Tx = new double[nz];
	const pointSet* mesh = M11.getMesh();
	double** data = M11.getData();

	int counter = 0;
	// Add Matrix 11
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = i;
			Tj[counter] = mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M12
	data = M12.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = i;
			Tj[counter] = n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M13
	data = M13.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = i;
			Tj[counter] = 2*n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M21
	data = M21.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = n+i;
			Tj[counter] = mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M22
	data = M22.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = n+i;
			Tj[counter] = n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M23
	data = M23.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = n+i;
			Tj[counter] = 2*n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M31
	data = M31.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = 2*n+i;
			Tj[counter] = mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M32
	data = M32.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = 2*n+i;
			Tj[counter] = n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	// Add Matrix M33
	data = M33.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
		{
			Ti[counter] = 2*n+i;
			Tj[counter] = 2*n+mesh->getIDX(i,j);
			Tx[counter] = data[i][j];
			counter++;
		}
	}

	if (counter >= nz )
		throw (10000);

	umfpack_di_triplet_to_col(this->n, this->n, counter, Ti, Tj, Tx, col, row, val, (int*) NULL);

	// Add RHS
	for (int i=0; i<n; i++)
	{
		rhs[i] = M11.getRhs(i) ;
		rhs[n+i] =  M22.getRhs(i);
		rhs[2*n+i] = M33.getRhs(i);

	}

	delete[] Ti;
	delete[] Tj;
	delete[] Tx;

}




void matrixUmfpack::solve(Array& Res)
{
	if ( Res.getSize() != n )
		throw "Error matrixUmfpack:: solve(Res): sizes differ!";
	double *null = (double*) NULL;
	void *Symbolic, *Numeric;
	double res[n];
	(void) umfpack_di_symbolic(n, n, col, row, val, &Symbolic, null, null);
	(void) umfpack_di_numeric(col, row, val, Symbolic, &Numeric, null, null);
	umfpack_di_free_symbolic(&Symbolic);
	(void) umfpack_di_solve(UMFPACK_A, col, row, val, res, rhs, Numeric, null, null);
	umfpack_di_free_numeric(&Numeric);

	for (int i=0; i<n; i++)
		Res(i) = res[i];
}

void matrixUmfpack::solve(Array& Res1, Array& Res2)
{
	if ( Res1.getSize() != n/2 )
		throw "Error matrixUmfpack:: solve(R1,R2): sizes differ!";
	double *null = (double*) NULL;
	void *Symbolic, *Numeric;
	double res[n];
	(void) umfpack_di_symbolic(n, n, col, row, val, &Symbolic, null, null);
	(void) umfpack_di_numeric(col, row, val, Symbolic, &Numeric, null, null);
	umfpack_di_free_symbolic(&Symbolic);
	(void) umfpack_di_solve(UMFPACK_A, col, row, val, res, rhs, Numeric, null, null);
	umfpack_di_free_numeric(&Numeric);

	for (int i=0; i<n/2; i++)
	{
		Res1(i) = res[i];
		Res2(i) = res[n/2+i];
	}
}


void matrixUmfpack::solve(Array& Res1, Array& Res2, Array& Res3)
{
	if ( Res1.getSize() != n/3 )
		throw "Error matrixUmfpack:: solve(R1,R2): sizes differ!";
	double *null = (double*) NULL;
	double info[UMFPACK_INFO];
	void *Symbolic, *Numeric;
	double res[n];
	(void) umfpack_di_symbolic(n, n, col, row, val, &Symbolic, null, null);
	(void) umfpack_di_numeric(col, row, val, Symbolic, &Numeric, null, info);
	umfpack_di_free_symbolic(&Symbolic);
	(void) umfpack_di_solve(UMFPACK_A, col, row, val, res, rhs, Numeric, null, info);
	umfpack_di_free_numeric(&Numeric);

	for (int i=0; i<n/3; i++)
	{
		Res1(i) = res[i];
		Res2(i) = res[n/3+i];
		Res3(i) = res[2*n/3+i];
	}
}


void matrixUmfpack::solve(scalarField& Res)
{
	Array& vec = Res();
	solve(vec);
}

void matrixUmfpack::solve(scalarField& Res1, scalarField& Res2)
{
	Array& vec1 = Res1();
	Array& vec2 = Res2();
	solve(vec1, vec2);
}

void matrixUmfpack::solve(scalarField& Res1, scalarField& Res2, scalarField& Res3)
{
	Array& vec1 = Res1();
	Array& vec2 = Res2();
	Array& vec3 = Res3();
	solve(vec1, vec2, vec3);
}


void matrixUmfpack::operator=(const matrixUmfpack& mat)
{
	if (n<0 && nz<0)
	{
		n  = mat.n;
		nz = mat.nz;
		col = new int[n+1];
		row = new int[nz];
		val = new double[nz];
		rhs = new double[n];
	}
	if (n!=mat.n | nz!=mat.nz)
		throw mfpmExcept(200);
	for (int i=0; i<n+1; i++)
		col[i] = mat.col[i];
	for (int i=0; i<n; i++)
		rhs[i] = mat.rhs[i];
	for (int i=0; i<nz; i++)
	{
		row[i] = mat.row[i];
		val[i] = mat.val[i];
	}

}



matrixUmfpack matrixUmfpack::transpose()
{
	matrixUmfpack mat(n, nz);
	for (int i=0; i<n; i++)
		mat.rhs[i] = rhs[i];

	int *null = (int*) NULL;
	(void) umfpack_di_transpose(n, n, col, row, val, null, null, mat.col, mat.row, mat.val);

	return mat;
}



void matrixUmfpack::calcDeriv
(
    const matrixUmfpack& matNew,
    const matrixUmfpack& matOld,
    double eps
)
{
	if (n!=matNew.n | nz!=matNew.nz
	        | n!=matOld.n | nz!=matOld.nz)
		throw mfpmExcept(200);
	for (int i=0; i<n+1; i++)
		col[i] = matNew.col[i];
	for (int i=0; i<n; i++)
		rhs[i] = matNew.rhs[i];
	for (int i=0; i<nz; i++)
	{
		row[i] = matNew.row[i];
		val[i] = (matNew.val[i] - matOld.val[i])/eps;
	}
}



double* matrixUmfpack::getRhs()
{
	return rhs;
}
