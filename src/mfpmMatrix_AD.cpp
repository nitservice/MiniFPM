#include "mfpmMatrix_AD.h"

mfpmMatrix_AD::mfpmMatrix_AD(const pointSet& meshIn):
	mesh( &meshIn )
{
	this->n    = meshIn.nbPoints();
	this->data = new adouble*[n];
	for (int i=0; i<n; i++)
		data[i] = new adouble[meshIn.getMaxPoints()];
	this->rhs  = new Array_AD(this->n);
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<meshIn.getMaxPoints(); j++)
			data[i][j] = 0.0;
		rhs->operator()(i) = 0.0;
	}

}


mfpmMatrix_AD::mfpmMatrix_AD(const mfpmMatrix_AD& mat):
	mesh( mat.getMesh() )
{
	this->n    = mesh->nbPoints();
	this->data = new adouble*[n];
	for (int i=0; i<n; i++)
		data[i] = new adouble[mesh->getMaxPoints()];
	this->rhs  = new Array_AD(this->n);
	this->operator=(mat);
}

mfpmMatrix_AD::mfpmMatrix_AD(const mfpmMatrix& mat):
	mesh( mat.getMesh() )
{
	this->n    = mesh->nbPoints();
	this->data = new adouble*[n];
	for (int i=0; i<n; i++)
		data[i] = new adouble[mesh->getMaxPoints()];
	this->rhs  = new Array_AD(this->n);
	this->operator=(mat);
}


mfpmMatrix_AD::mfpmMatrix_AD()
{
	data = NULL;
	rhs  = NULL;
	n    = 0;
}

mfpmMatrix_AD::~mfpmMatrix_AD()
{
	if (n>0)
	{
		for (int i=0; i<n; i++)
			delete[] data[i];
		delete[] data;
	}
	delete rhs;
}


void mfpmMatrix_AD::operator=(const mfpmMatrix_AD& mat)
{
	if (n<=0)
	{
		cout << "Allocate matrix " << endl;
		this->mesh = mat.getMesh();
		this->n    = mesh->nbPoints();
		this->data = new adouble*[n];
		for (int i=0; i<n; i++)
			data[i] = new adouble[mesh->getMaxPoints()];
		this->rhs  = new Array_AD(this->n);
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
		this->data = new adouble*[n];
		for (int i=0; i<n; i++)
			data[i] = new adouble[mesh->getMaxPoints()];
		this->rhs  = new Array_AD(this->n);
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


void mfpmMatrix_AD::operator=(const mfpmMatrix& mat)
{
	if (n<=0)
	{
		cout << "Allocate matrix " << endl;
		this->mesh = mat.getMesh();
		this->n    = mesh->nbPoints();
		this->data = new adouble*[n];
		for (int i=0; i<n; i++)
			data[i] = new adouble[mesh->getMaxPoints()];
		this->rhs  = new Array_AD(this->n);
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
		this->data = new adouble*[n];
		for (int i=0; i<n; i++)
			data[i] = new adouble[mesh->getMaxPoints()];
		this->rhs  = new Array_AD(this->n);
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





adouble mfpmMatrix_AD::operator()(int row, int col) const
{
	int idn=0;
	for (; idn<mesh->getNIDX(row); idn++)
		if (mesh->getIDX(row,idn) == col)
			break;
	if (idn<mesh->getNIDX(row))
		return data[row][idn];
	else
		return 0.0;
}


void mfpmMatrix_AD::setRow(int row, Array_AD val)
{
	int m = mesh->getNIDX(row);
	if (val.getSize() != m)
		throw "mfpmMatrix_AD error:: setRow: invalid sizes!";
	for (int i=0; i<m; i++)
		data[row][i] = val(i);
}


void mfpmMatrix_AD::setDiag(int row, adouble& val)
{
	for (int j=0; j<mesh->getNIDX(row); j++)
		if (mesh->getIDX(row,j) == row)
			data[row][j] = val;
		else
			data[row][j] = 0.0;
}

void mfpmMatrix_AD::addRow(int row, Array_AD val)
{
	int m = mesh->getNIDX(row);
	if (val.getSize() != m)
		throw "mfpmMatrix_AD error:: setRow: invalid sizes!";
	for (int i=0; i<m; i++)
		data[row][i] = data[row][i] + val(i);
}

Array_AD mfpmMatrix_AD::getRow(int row) const
{
	int m = mesh->getNIDX(row);
	Array_AD res(m);
	for (int i=0; i<m; i++)
	{
		res(i) = data[row][i];
	}

	return res;
}


void mfpmMatrix_AD::clearBnd()
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




void mfpmMatrix_AD::setDirichlet(int bnd, double val)
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


void mfpmMatrix_AD::setDirichlet(int bnd, const Array_AD& val)
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


void mfpmMatrix_AD::setNeumann(int bnd, double val)
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

void mfpmMatrix_AD::setNeumann(int bnd, const Array_AD& val)
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


void mfpmMatrix_AD::setRobin(int bnd, double eps, double val)
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


void mfpmMatrix_AD::setRobin(int bnd, double eps, const Array_AD& val)
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



Array_AD mfpmMatrix_AD::operator*(const Array_AD& x) const
{
	if (x.getSize() != n)
		throw "mfpmMatrix_AD error:: operator*: invalid sizse!";
	Array_AD Result(n);
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


Array_AD mfpmMatrix_AD::operator*(const scalarField& x) const
{
	if (x.getSize() != n)
		throw "mfpmMatrix_AD error:: operator*: invalid sizse!";
	Array_AD Result(n);
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


mfpmMatrix_AD mfpmMatrix_AD::operator*(const adouble& d) const
{
	mfpmMatrix_AD mat(*this);
	Array_AD* rhs = mat.getRhs();
	adouble** md = mat.getData();
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



mfpmMatrix_AD operator*(const adouble& d, const mfpmMatrix_AD& mat)
{
	return mat*d;
}




mfpmMatrix_AD mfpmMatrix_AD::operator+(const Array_AD& V) const
{
	if (V.getSize() != n)
		throw "Error mfpmMatrix_AD:: operator+: wrong dimensions";
	mfpmMatrix_AD mat = *this;
	Array_AD* mr = mat.getRhs();
	for (int i=0; i<n; i++)
		mr->operator()(i) = rhs->value(i) - V.value(i);
	return mat;
}



mfpmMatrix_AD mfpmMatrix_AD::operator+(const mfpmMatrix_AD& mat) const
{
	if (mat.getMesh() != mesh)
		throw "Error mfpmMatrix_AD:: operator+: wrong dimensions";
	mfpmMatrix_AD tmp = *this;
	Array_AD*  tmr = tmp.getRhs();
	adouble** tmd = tmp.getData();
	Array_AD*  imr = mat.getRhs();
	adouble** imd = mat.getData();
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<mesh->getNIDX(i); j++)
			tmd[i][j] = data[i][j] + imd[i][j];
		tmr->operator()(i) = rhs->value(i) + imr->value(i);
	}
	return tmp;
}


mfpmMatrix_AD mfpmMatrix_AD::operator+(const adouble& d) const
{
	mfpmMatrix_AD mat = *this;
	Array_AD* mr = mat.getRhs();
	for (int i=0; i<n; i++)
		mr->operator()(i) = rhs->value(i) - d;
	return mat;
}


mfpmMatrix_AD operator*(const Array_AD& d, const mfpmMatrix_AD& mat)
{
	mfpmMatrix_AD rmat(*(mat.getMesh()));
	Array_AD* mr = rmat.getRhs();
	for (int i=0; i<mat.getSize(); i++)
	{
		Array_AD row = mat.getRow(i);
		row = d(i)*row;
		rmat.setRow(i, row);
		(*mr)(i) = mat.getRhs(i)* d(i);
	}
	return rmat;
}


mfpmMatrix_AD operator+(const adouble& d, const mfpmMatrix_AD& mat)
{
	return mat + d;
}

mfpmMatrix_AD operator+(const Array_AD& V, const mfpmMatrix_AD& mat)
{
	return mat + V;
}


mfpmMatrix_AD mfpmMatrix_AD::operator-(const adouble& d) const
{
	return (*this + -d);
}

mfpmMatrix_AD mfpmMatrix_AD::operator-(const Array_AD& V) const
{
	return (*this + -1.0*V);
}

mfpmMatrix_AD mfpmMatrix_AD::operator-(const mfpmMatrix_AD& mat) const
{
	adouble tmp = 1.0;
	return (*this + -tmp*mat);
}


mfpmMatrix_AD operator-(const adouble& d, const mfpmMatrix_AD& mat)
{
	adouble tmp = 1.0;
	return d + -tmp*mat;
}

mfpmMatrix_AD operator-(const Array_AD& V, const mfpmMatrix_AD& mat)
{
	adouble tmp = 1.0;
	return V + -tmp*mat;
}


mfpmMatrix_AD operator==(const mfpmMatrix_AD& mat, const Array_AD& V)
{
	return mat - V;
}

mfpmMatrix_AD operator==(const mfpmMatrix_AD& mat, const adouble& d)
{
	return mat - d;
}





void mfpmMatrix_AD::regularise(double factor)
{
	for (int i=0; i<n; i++)
		//if (mesh->getBND(i)==0)
		for (int j=0; j<mesh->getNIDX(i); j++)
			data[i][j] = data[i][j] + factor;

}


void mfpmMatrix_AD::normalise()
{
	adouble diagVal = 0;
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
			throw "Error mfpmMatrix_AD:: normalise():: zero diagVal!";
		for (int j=0; j<mesh->getNIDX(i); j++)
			data[i][j] = data[i][j] / diagVal;
		(*rhs)(i) = (*rhs)(i) / diagVal;
	}
}



void mfpmMatrix_AD::addRhs(int i, adouble& val)
{
	(*rhs)(i) += val;
}

void mfpmMatrix_AD::setRhs(int i, adouble& val)
{
	(*rhs)(i) = val;
}
