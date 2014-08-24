#include "mfpmMatrix.h"
#include "pointSet.h"
#include "scalarField.h"
#include "primitive_operators.h"

#include "matrixUmfpack.h"

mfpmMatrix::mfpmMatrix(const pointSet& meshIn):
	mesh( &meshIn )
{
	this->n    = meshIn.nbPoints();
	this->data = new double*[n];
	int maxPoints = meshIn.getMaxPoints();
	for (int i=0; i<n; i++)
		data[i] = new double[meshIn.getNIDX(i)];
	this->rhs  = new Array(this->n);
//#pragma omp parallel for
	for (int i=0; i<n; i++)
	{
//		for (int j=0; j<meshIn.getMaxPoints(); j++)
		for (int j=0; j<meshIn.getNIDX(i); j++)
			data[i][j] = 0.0;
		rhs->operator()(i) = 0.0;
	}

}


mfpmMatrix::mfpmMatrix(const mfpmMatrix& mat):
	mesh( mat.getMesh() )
{
	this->n    = mesh->nbPoints();
	this->data = new double*[n];
	for (int i=0; i<n; i++)
		data[i] = new double[mesh->getMaxPoints()];
	this->rhs  = new Array(this->n);
	this->operator=(mat);
}


mfpmMatrix::mfpmMatrix()
{
	data = NULL;
	rhs  = NULL;
	n    = 0;
}

mfpmMatrix::~mfpmMatrix()
{
	if (n>0)
	{
		for (int i=0; i<n; i++)
			delete[] data[i];
		delete[] data;
	}
	delete rhs;
}


void mfpmMatrix::operator=(const mfpmMatrix& mat)
{
	if (n<=0)
	{
		cout << "Allocate matrix " << endl;
		this->mesh = mat.getMesh();
		this->n    = mesh->nbPoints();
		this->data = new double*[n];
		for (int i=0; i<n; i++)
			data[i] = new double[mesh->getMaxPoints()];
		this->rhs  = new Array(this->n);
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<mesh->getNIDX(i); j++)
				data[i][j] = mat.getData(i,j);
			rhs->operator()(i) = mat.getRhs(i);
		}
	}
	else if (this->n != mat.getSize())
	{
		cout << "Reshape matrix!" << endl;
		if (data != NULL)
		{
			for (int i=0; i<mesh->nbPoints(); i++)
				delete[] data[i];
			delete[] data;
		}
		delete rhs;
		this->mesh = mat.getMesh();
		this->n    = mesh->nbPoints();
		this->data = new double*[n];
		for (int i=0; i<n; i++)
			data[i] = new double[mesh->getMaxPoints()];
		this->rhs  = new Array(this->n);
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<mesh->getNIDX(i); j++)
				data[i][j] = mat.getData(i,j);
			rhs->operator()(i) = mat.getRhs(i);
		}
	}
	else
	{
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<mesh->getNIDX(i); j++)
				data[i][j] = mat.getData(i,j);
			rhs->operator()(i) = mat.getRhs(i);
		}
	}

}





double mfpmMatrix::operator()(int row, int col) const
{
	int idn=0;
	bool found = false;
	for (; idn<mesh->getNIDX(row); idn++)
		if (mesh->getIDX(row,idn) == col)
		{
			found = true;
			break;
		}
	if (found)
		return data[row][idn];
	else
		return 0.0;
}


void mfpmMatrix::setRow(int row, Array val)
{
	int m = mesh->getNIDX(row);
	if (val.getSize() != m)
		throw "mfpmMatrix error:: setRow: invalid sizes!";
	for (int i=0; i<m; i++)
		data[row][i] = val(i);
}


void mfpmMatrix::setDiag(int row, double val)
{
	for (int j=0; j<mesh->getNIDX(row); j++)
		if (mesh->getIDX(row,j) == row)
			data[row][j] = val;
		else
			data[row][j] = 0.0;
}

void mfpmMatrix::addRow(int row, Array val)
{
	int m = mesh->getNIDX(row);
	if (val.getSize() != m)
		throw "mfpmMatrix error:: setRow: invalid sizes!";
	for (int i=0; i<m; i++)
		data[row][i] = data[row][i] + val(i);
}

Array mfpmMatrix::getRow(int row) const
{
	int m = mesh->getNIDX(row);
	Array res(m);
	for (int i=0; i<m; i++)
	{
		res(i) = data[row][i];
	}

	return res;
}


void mfpmMatrix::clearBnd()
{
	for (int i=0; i<n; i++)
	{
		if (mesh->getBND(i) != 0)
		{
			for (int j=0; j<mesh->getNIDX(i); j++)
				data[i][j] = 0.0;
			(*rhs)(i) = 0.0;
		}
	}
}




void mfpmMatrix::setDirichlet(int bnd, double val)
{
	for (int i=0; i<n; i++)
		if (mesh->getBND(i) == bnd)
		{
			for (int j=0; j<mesh->getNIDX(i); j++)
				if (mesh->getIDX(i,j) == i)
					data[i][j] = 1.0;
				else
					data[i][j] = 0.0;
			rhs->operator()(i) = val;
		}
}


void mfpmMatrix::setDirichlet(int bnd, const Array& val)
{
	for (int i=0; i<n; i++)
		if (mesh->getBND(i) == bnd)
		{
			for (int j=0; j<mesh->getNIDX(i); j++)
				if (mesh->getIDX(i,j) == i)
					data[i][j] = 1.0;
				else
					data[i][j] = 0.0;
			rhs->operator()(i) = val(i);
		}
}


void mfpmMatrix::setNeumann(int bnd, double val)
{
	for (int i=0; i<n; i++)
		if (mesh->getBND(i) == bnd)
		{
			int nidx = mesh->getNIDX(i);
			Array tmp(nidx);
			int info = neumannMatEntry(*mesh, i, tmp);
			if (info < 0)
			{
				info = neumannMatEntry(*mesh, i, tmp, true);
				if (info < 0)
					cout << "Error in Neuman condition!!!!!!!!!!!!!!" << endl;
			}

			for (int j=0; j<nidx; j++)
				data[i][j] = tmp(j);
			rhs->operator()(i) = val;
		}
}

void mfpmMatrix::setNeumann(int bnd, const Array& val)
{
	for (int i=0; i<n; i++)
		if (mesh->getBND(i) == bnd)
		{
			int nidx = mesh->getNIDX(i);
			Array tmp(nidx);
			int info = neumannMatEntry(*mesh, i, tmp);
			if (info < 0)
			{
				info = neumannMatEntry(*mesh, i, tmp, true);
				if (info < 0)
					cout << "Error in Neuman condition!!!!!!!!!!!!!!" << endl;
			}
			for (int j=0; j<nidx; j++)
				data[i][j] = tmp(j);
			rhs->operator()(i) = val(i);
		}
}


void mfpmMatrix::setRobin(int bnd, double eps, double val)
{
	// Robin boundary condition, i.e.
	//   eps*grad(u)*n + u = val   on bnd
	setDirichlet(bnd, val);
	for (int i=0; i<n; i++)
		if (mesh->getBND(i) == bnd)
		{
			int nidx = mesh->getNIDX(i);
			Array tmp(nidx);
			int info = neumannMatEntry(*mesh, i, tmp);
			if (info < 0)
			{
				info = neumannMatEntry(*mesh, i, tmp, true);
				if (info < 0)
					cout << "Error in Neuman condition!!!!!!!!!!!!!!" << endl;
			}
			for (int j=0; j<nidx; j++)
				data[i][j] += eps*tmp(j);
		}
}


void mfpmMatrix::setRobin(int bnd, double eps, const Array& val)
{
	setDirichlet(bnd, val);
	for (int i=0; i<n; i++)
		if (mesh->getBND(i) == bnd)
		{
			int nidx = mesh->getNIDX(i);
			Array tmp(nidx);
			int info = neumannMatEntry(*mesh, i, tmp);
			if (info < 0)
			{
				info = neumannMatEntry(*mesh, i, tmp, true);
				if (info < 0)
					cout << "Error in Neuman condition!!!!!!!!!!!!!!" << endl;
			}
			for (int j=0; j<nidx; j++)
				data[i][j] += eps*tmp(j);
		}
}



Array mfpmMatrix::operator*(const Array& x) const
{
	if (x.getSize() != n)
		throw "mfpmMatrix error:: operator*: invalid sizse!";
	Array Result(n);
	for (int row=0; row<n; row++)
	{
		Result(row) = 0.0;
		for (int col=0; col<mesh->getNIDX(row); col++)
		{
			Result(row) += x.value(mesh->getIDX(row,col))*data[row][col];
		}
	}
	return Result;
}


Array mfpmMatrix::operator*(const scalarField& x) const
{
	if (x.getSize() != n)
		throw "mfpmMatrix error:: operator*: invalid sizse!";
	Array Result(n);
	for (int row=0; row<n; row++)
	{
		Result(row) = 0.0;
		for (int col=0; col<mesh->getNIDX(row); col++)
		{
			Result(row) += x.value(mesh->getIDX(row,col))*data[row][col];
		}
	}
	return Result;
}


mfpmMatrix mfpmMatrix::operator*(const double d) const
{
	mfpmMatrix mat(*this);
	Array* rhs = mat.getRhs();
	double** md = mat.getData();
	for (int row=0; row<n; row++)
	{
		for (int col=0; col<mesh->getNIDX(row); col++)
		{
			md[row][col] = data[row][col] * d;
		}
		rhs->operator()(row)   = (this->rhs)->operator()(row) * d;

	}
	return mat;
}



mfpmMatrix operator*(const double d, const mfpmMatrix& mat)
{
	return mat*d;
}




mfpmMatrix mfpmMatrix::operator+(const Array& V) const
{
	if (V.getSize() != n)
		throw "Error mfpmMatrix:: operator+: wrong dimensions";
	mfpmMatrix mat = *this;
	Array* mr = mat.getRhs();
	for (int i=0; i<n; i++)
		mr->operator()(i) = rhs->value(i) - V.value(i);
	return mat;
}



mfpmMatrix mfpmMatrix::operator+(const mfpmMatrix& mat) const
{
	if (mat.getMesh() != mesh)
		throw "Error mfpmMatrix:: operator+: wrong dimensions";
	mfpmMatrix tmp = *this;
	Array*  tmr = tmp.getRhs();
	double** tmd = tmp.getData();
	Array*  imr = mat.getRhs();
	double** imd = mat.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
			tmd[i][j] = data[i][j] + imd[i][j];
		tmr->operator()(i) = rhs->value(i) + imr->value(i);
	}
	return tmp;
}


mfpmMatrix mfpmMatrix::operator+(const double d) const
{
	mfpmMatrix mat = *this;
	Array* mr = mat.getRhs();
	for (int i=0; i<n; i++)
		mr->operator()(i) = rhs->value(i) - d;
	return mat;
}


mfpmMatrix operator*(const Array& d, const mfpmMatrix& mat)
{
	mfpmMatrix rmat(*(mat.getMesh()));
	Array* mr = rmat.getRhs();
	for (int i=0; i<mat.getSize(); i++)
	{
		Array row = mat.getRow(i);
		row = d(i)*row;
		rmat.setRow(i, row);
		(*mr)(i) = mat.getRhs(i)* d(i);
	}
	return rmat;
}


mfpmMatrix operator+(const double d, const mfpmMatrix& mat)
{
	return mat + d;
}

mfpmMatrix operator+(const Array& V, const mfpmMatrix& mat)
{
	return mat + V;
}


mfpmMatrix mfpmMatrix::operator-(const double d) const
{
	return (*this + -d);
}

mfpmMatrix mfpmMatrix::operator-(const Array& V) const
{
	return (*this + -1.0*V);
}

mfpmMatrix mfpmMatrix::operator-(const mfpmMatrix& mat) const
{
	return (*this + -1.0*mat);
}


mfpmMatrix operator-(const double d, const mfpmMatrix& mat)
{
	return d + -1.0*mat;
}

mfpmMatrix operator-(const Array& V, const mfpmMatrix& mat)
{
	return V + -1.0*mat;
}


mfpmMatrix operator==(const mfpmMatrix& mat, const Array& V)
{
	return mat - V;
}

mfpmMatrix operator==(const mfpmMatrix& mat, const double d)
{
	return mat - d;
}


void mfpmMatrix::solveUMFPACK(scalarField& Res)
{
	matrixUmfpack mat(n, n*mesh->getMaxPoints());
	mat.setMatrix(*this);
	mat.solve(Res);
}



void mfpmMatrix::regularise(double factor)
{
	for (int i=0; i<n; i++)
		//if (mesh->getBND(i)==0)
		for (int j=0; j<mesh->getNIDX(i); j++)
			data[i][j] = data[i][j] + factor;

}


void mfpmMatrix::normalise()
{
	double diagVal = 0;
	for (int i=0; i<n; i++)
	{
		diagVal = 0.0;
		for (int j=0; j<mesh->getNIDX(i); j++)
			if (mesh->getIDX(i,j) == i)
			{
				diagVal = data[i][j];
				break;
			}
		if (diagVal == 0.0)
			throw "Error mfpmMatrix:: normalise():: zero diagVal!";
		for (int j=0; j<mesh->getNIDX(i); j++)
			data[i][j] = data[i][j] / diagVal;
		(*rhs)(i) = (*rhs)(i) / diagVal;
	}
}



void mfpmMatrix::addRhs(int i, double val)
{
	(*rhs)(i) += val;
}

void mfpmMatrix::setRhs(int i, double val)
{
	(*rhs)(i) = val;
}
