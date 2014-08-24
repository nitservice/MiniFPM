#include "matrixIML.h"

#include <comprow_double.h>
#include <diagpre_double.h>
#include <icpre_double.h>
//#include <ilupre_double.h>
#include <mvblasd.h>

#include <bicgstab.h>


int solve3x3(
    const mfpmMatrix& a11,
    const mfpmMatrix& a12,
    const mfpmMatrix& a13,
    const mfpmMatrix& a21,
    const mfpmMatrix& a22,
    const mfpmMatrix& a23,
    const mfpmMatrix& a31,
    const mfpmMatrix& a32,
    const mfpmMatrix& a33,
    Array& res1,
    Array& res2,
    Array& res3
)
{

	int n = a11.getSize();
	const pointSet* mesh = a11.getMesh();
	// Fill matrix
	int NNZ = 9*n*mesh->getMaxPoints();
	double* tmpVal = new double[NNZ];
	int*    tmpPos = new int[NNZ];
	int*    tmpRow = new int[3*n+1];
	int cnt = 0;

	for (int i=0; i<NNZ; i++)
	{
		tmpPos[i] = 0;
		tmpVal[i] = 0;
	}
	for (int i=0; i<3*n+1; i++)
	{
		tmpRow[i] = 0;
	}

	for (int i=0; i<n; i++)
	{
		tmpRow[i] = cnt;
		Array row = a11.getRow(i);
		for (int j=0; j<row.getSize(); j++)
		{
			tmpVal[cnt] = row(j);
			tmpPos[cnt] = mesh->getIDX(i,j);
			cnt++;
		}
		Array row2 = a12.getRow(i);
		for (int j=0; j<row2.getSize(); j++)
		{
			tmpVal[cnt] = row2(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + n;
			cnt++;
		}
		Array row3 = a13.getRow(i);
		for (int j=0; j<row3.getSize(); j++)
		{
			tmpVal[cnt] = row3(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + 2*n;
			cnt++;
		}
	}

	for (int i=0; i<n; i++)
	{
		tmpRow[i+n] = cnt;
		Array row = a21.getRow(i);
		for (int j=0; j<row.getSize(); j++)
		{
			tmpVal[cnt] = row(j);
			tmpPos[cnt] = mesh->getIDX(i,j);
			cnt++;
		}
		Array row2 = a22.getRow(i);
		for (int j=0; j<row2.getSize(); j++)
		{
			tmpVal[cnt] = row2(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + n;
			cnt++;
		}
		Array row3 = a23.getRow(i);
		for (int j=0; j<row3.getSize(); j++)
		{
			tmpVal[cnt] = row3(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + 2*n;
			cnt++;
		}
		int m = i + n;
	}

	for (int i=0; i<n; i++)
	{
		tmpRow[i+2*n] = cnt;
		Array row = a31.getRow(i);
		for (int j=0; j<row.getSize(); j++)
		{
			tmpVal[cnt] = row(j);
			tmpPos[cnt] = mesh->getIDX(i,j);
			cnt++;
		}
		Array row2 = a32.getRow(i);
		for (int j=0; j<row2.getSize(); j++)
		{
			tmpVal[cnt] = row2(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + n;
			cnt++;
		}
		Array row3 = a33.getRow(i);
		for (int j=0; j<row3.getSize(); j++)
		{
			tmpVal[cnt] = row3(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + 2*n;
			cnt++;
		}
		int m = i + 2*n;
	}
	tmpRow[3*n] = cnt;

	CompRow_Mat_double SysMat(3*n, 3*n, cnt	, tmpVal, tmpRow, tmpPos);
	DiagPreconditioner_double diagPre(SysMat);
	//CompRow_ILUPreconditioner_double diagPre(SysMat);
	VECTOR_double res(3*n);
	VECTOR_double rhs(3*n);
	for (int i=0; i<n; i++)
	{
		rhs(i) = a11.getRhs(i);
		rhs(i+n) = a22.getRhs(i);
		rhs(i+2*n) = a33.getRhs(i);
		res(i) = 0.0;
		res(i+n) = 0.0;
		res(i+2*n) = 0.0;
	}

	int maxit = 2000;
	double tol = 1e-5;
	int result = BiCGSTAB(SysMat, res, rhs, diagPre, maxit, tol);
	if (result != 0)
	{
		cout << "Status:" << endl;
		cout << " -> flag: " << result << endl;
		cout << " -> It:   " << maxit << endl;
		cout << " -> tol:  " << tol << endl;
	}
	for (int i=0; i<n; i++)
	{
		res1(i) = res(i);
		res2(i) = res(i+n);
		res3(i) = res(i+2*n);
	}

	delete[] tmpPos;
	delete[] tmpVal;
	delete[] tmpRow;

	return result;
}



int solve2x2(
    const mfpmMatrix& a11,
    const mfpmMatrix& a12,
    const mfpmMatrix& a21,
    const mfpmMatrix& a22,
    Array& res1,
    Array& res2
)
{

	int n = a11.getSize();
	const pointSet* mesh = a11.getMesh();
	// Fill matrix
	int NNZ = 4*n*mesh->getMaxPoints();
	double* tmpVal = new double[NNZ];
	int*    tmpPos = new int[NNZ];
	int*    tmpRow = new int[2*n+1];
	int cnt = 0;

	for (int i=0; i<NNZ; i++)
	{
		tmpPos[i] = 0;
		tmpVal[i] = 0;
	}
	for (int i=0; i<2*n+1; i++)
	{
		tmpRow[i] = 0;
	}

	for (int i=0; i<n; i++)
	{
		tmpRow[i] = cnt;
		Array row = a11.getRow(i);
		for (int j=0; j<row.getSize(); j++)
		{
			tmpVal[cnt] = row(j);
			tmpPos[cnt] = mesh->getIDX(i,j);
			cnt++;
		}
		Array row2 = a12.getRow(i);
		for (int j=0; j<row2.getSize(); j++)
		{
			tmpVal[cnt] = row2(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + n;
			cnt++;
		}
	}

	for (int i=0; i<n; i++)
	{
		tmpRow[i+n] = cnt;
		Array row = a21.getRow(i);
		for (int j=0; j<row.getSize(); j++)
		{
			tmpVal[cnt] = row(j);
			tmpPos[cnt] = mesh->getIDX(i,j);
			cnt++;
		}
		Array row2 = a22.getRow(i);
		for (int j=0; j<row2.getSize(); j++)
		{
			tmpVal[cnt] = row2(j);
			tmpPos[cnt] = mesh->getIDX(i,j) + n;
			cnt++;
		}
		int m = i + n;
	}
	tmpRow[2*n] = cnt;

	CompRow_Mat_double SysMat(2*n, 2*n, cnt	, tmpVal, tmpRow, tmpPos);
	DiagPreconditioner_double diagPre(SysMat);
	//CompRow_ILUPreconditioner_double diagPre(SysMat);
	VECTOR_double res(2*n);
	VECTOR_double rhs(2*n);
	for (int i=0; i<n; i++)
	{
		rhs(i) = a11.getRhs(i);
		rhs(i+n) = a22.getRhs(i);
		res(i) = 0.0;
		res(i+n) = 0.0;
	}

	int maxit = 2000;
	double tol = 1e-5;
	int result = BiCGSTAB(SysMat, res, rhs, diagPre, maxit, tol);
	if (result != 0)
	{
		cout << "Status:" << endl;
		cout << " -> flag: " << result << endl;
		cout << " -> It:   " << maxit << endl;
		cout << " -> tol:  " << tol << endl;
	}
	for (int i=0; i<n; i++)
	{
		res1(i) = res(i);
		res2(i) = res(i+n);
	}

	delete[] tmpPos;
	delete[] tmpVal;
	delete[] tmpRow;

	return result;
}




int solve1x1(
    const mfpmMatrix& A,
    Array& res1
)
{

	int n = A.getSize();
	const pointSet* mesh = A.getMesh();
	// Fill matrix
	int NNZ = n*mesh->getMaxPoints();
	double* tmpVal = new double[NNZ];
	int*    tmpPos = new int[NNZ];
	int*    tmpRow = new int[n+1];
	int cnt = 0;

	for (int i=0; i<n; i++)
	{
		tmpRow[i] = cnt;
		Array row = A.getRow(i);
		for (int j=0; j<row.getSize(); j++)
		{
			tmpVal[cnt] = row(j);
			tmpPos[cnt] = mesh->getIDX(i,j);
			cnt++;
		}
	}

	tmpRow[n] = cnt;

	CompRow_Mat_double SysMat(n, n, cnt	, tmpVal, tmpRow, tmpPos);
	DiagPreconditioner_double diagPre(SysMat);
	VECTOR_double res(n);
	VECTOR_double rhs(n);
	for (int i=0; i<n; i++)
	{
		rhs(i) = A.getRhs(i);
		res(i) = 0.0;
	}

	int maxit = 2000;
	double tol = 1e-5;
	int result = BiCGSTAB(SysMat, res, rhs, diagPre, maxit, tol);
	if (result != 0)
	{
		cout << "Status:" << endl;
		cout << " -> flag: " << result << endl;
		cout << " -> It:   " << maxit << endl;
		cout << " -> tol:  " << tol << endl;
	}
	for (int i=0; i<n; i++)
	{
		res1(i) = res(i);
	}

	delete[] tmpPos;
	delete[] tmpVal;
	delete[] tmpRow;

	return result;
}
