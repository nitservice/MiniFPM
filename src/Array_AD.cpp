#include "Array_AD.h"
#include <iostream>
#include <cmath>
#include <adolc/adolc.h>
#include "Array.h"


Array_AD::Array_AD(int n)
{
	size = n;
	data.resize(n,0.0);
}

Array_AD::Array_AD(const Array_AD& V)
{
	size = V.getSize();
	data.resize(size);
	for ( int i=0; i<size; i++)
		data[i] = V.value(i);
}

Array_AD::Array_AD(const Array& V)
{
	size = V.getSize();
	data.resize(size);
	for ( int i=0; i<size; i++)
		data[i]=V.value(i);
}

Array_AD::Array_AD(int n, adouble* val)
{
	size = n;
	data.resize(size);
	for ( int i=0; i<size; i++)
		data[i]= val[i];
}

Array_AD::Array_AD()
{
	size = 0;
}

Array_AD::~Array_AD()
{
	data.clear();
}


adouble& Array_AD::operator()(int pos)
{
	return data[pos];
}

adouble Array_AD::operator()(int pos) const
{
	return data[pos];
}



void Array_AD::operator=(const Array_AD& V)
{
	if (data.empty())
	{
		size = V.getSize();
		data.resize(size);
		for ( int i=0; i<size; i++)
			data[i] = V.value(i);
	}
	else
	{
		if (size!=V.getSize())
		{
			size = V.getSize();
			data.resize(size);
		}
		for (int i=0; i<size; i++)
			data[i] = V.value(i);
	}
}

void Array_AD::operator=(const Array& V)
{
	if (data.empty())
	{
		size = V.getSize();
		data.resize(size);
		for ( int i=0; i<size; i++)
			data[i] = V.value(i);
	}
	else
	{
		if (size!=V.getSize())
		{
			size = V.getSize();
			data.resize(size);
		}
		for (int i=0; i<size; i++)
			data[i] = V.value(i);
	}
}


void Array_AD::operator=(const adouble& val)
{
	if (data.empty())
	{
		throw "Error Array_AD:: operator=: object not initialised!";
	}
	else
	{
		for (int i=0; i<size; i++)
			data[i] = val;
	}
}


Array_AD Array_AD::operator+(const Array_AD& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array_AD class:: operator+: sizes differ";
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] + V.value(i);
	return A;
}


Array_AD Array_AD::operator+(const adouble& v) const
{
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] + v;
	return A;
}

Array_AD operator+(const adouble d, const Array_AD& V)
{
	return V+d;
}


Array_AD Array_AD::operator-(const Array_AD& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array_AD class:: operator+: sizes differ";
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] - V.value(i);
	return A;
}

Array_AD Array_AD::operator-(const adouble& v) const
{
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] - v;
	return A;
}

Array_AD operator-(const adouble& d, const Array_AD& V)
{
	Array_AD A(V.getSize());
	for (int i=0; i<V.getSize(); i++)
		A(i) = d - V.value(i);
	return A;
}


Array_AD operator-(const Array_AD& V)
{
	Array_AD A(V.getSize());
	for (int i=0; i<V.getSize(); i++)
		A(i) = - V.value(i);
	return A;
}



Array_AD Array_AD::operator*(const Array_AD& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array_AD class:: operator*: sizes differ";
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] * V.value(i);
	return A;
}


Array_AD Array_AD::operator*(const adouble& V) const
{
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] * V;
	return A;
}

Array_AD operator*(const adouble& d, const Array_AD& V)
{
	return V*d;
}



Array_AD Array_AD::operator/(const Array_AD& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array_AD class:: operator/: sizes differ";
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] / V.value(i);
	return A;
}

Array_AD Array_AD::operator/(const adouble& d) const
{
	Array_AD A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] / d;
	return A;
}


Array_AD operator/(const adouble& d, const Array_AD& V)
{
	Array_AD A(V.getSize());
	for (int i=0; i<V.getSize(); i++)
		A(i) = d / V.value(i);
	return A;
}

Array_AD Array_AD::operator<(const Array_AD& A) const
{
	Array_AD tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] < A(i) )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}


Array_AD Array_AD::operator<(adouble& d) const
{
	Array_AD tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] < d )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}


Array_AD Array_AD::operator>(const Array_AD& A) const
{
	Array_AD tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] > A(i) )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}


Array_AD Array_AD::operator>(adouble& d) const
{
	Array_AD tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] > d )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}



//   ADOLC operators

void Array_AD::operator<<=(const Array& V)
{
	if (data.empty())
	{
		size = V.getSize();
		data.resize(size);
		for ( int i=0; i<size; i++) {
			data[i] <<= V.value(i);
		}
	}
	else
	{
		if (size!=V.getSize())
		{
			size = V.getSize();
			data.resize(size);
		}
		for ( int i=0; i<size; i++) {
			adouble val;
			data[i] <<= V.value(i);
		}
	}
}


void Array_AD::operator>>=(Array& V)
{
	for (int i=0; i<size; i++)
	{
		data[i] >>= V(i);
	}
}


void Array_AD::setSize(int n)
{
	size = n;
	data.resize(n, 0.0);
}


void Array_AD::remove(int pn)
{
	if (pn<size)
	{
		vector<adouble>::iterator it;
		it = data.begin() + pn;
		data.erase(it);
		size=data.size();
	}
	else
		throw "Error Array_AD::remove: invalid point number";
}



void Array_AD::add(int pn, adouble& val)
{
	if (pn<=size & pn>=0)
	{
		data.insert(data.begin() + pn, val);
		size=data.size();
	}
	else
		throw "Error Array_AD::add(int, adouble): invalid point number";
}


void Array_AD::add(int pn)
{
	if (pn<=size & pn>=0)
	{
		data.insert(data.begin() + pn, 0.0);
		size=data.size();
	}
	else
		throw "Error Array_AD::add(int, adouble): invalid point number";
}


void Array_AD::append(adouble& val)
{
	data.push_back(val);
	size=data.size();
}


adouble minval(const Array_AD& vec)
{
	adouble mv = vec(0);
	for (int i=1; i<vec.size; i++)
		if (mv>vec(i))
			mv = vec(i);
	return mv;
}

adouble maxval(const Array_AD& vec)
{
	adouble mv = vec(0);
	for (int i=1; i<vec.size; i++)
		if (mv<vec(i))
			mv = vec(i);
	return mv;
}

void Array_AD::print()
{
	int n = getSize();
	cout << "Array_AD values: " << n << endl;
	cout << "============================================" << endl;
	for (int i=0; i<n; i++)
		cout << data[i] << endl;
	cout << "============================================" << endl;
}



/********************************
 * Math operators, e.g. exp, pow
 ********************************/
Array_AD exp(const Array_AD& A)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = exp(A(i));
	return tmp;
}

Array_AD sqrt(const Array_AD& A)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = sqrt(A(i));
	return tmp;
}

Array_AD sin(const Array_AD& A)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = sin(A(i));
	return tmp;
}

Array_AD cos(const Array_AD& A)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = cos(A(i));
	return tmp;
}

Array_AD pow(const Array_AD& A, int p)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = pow(A(i), p);
	return tmp;
}



Array_AD min(const Array_AD& A, const Array_AD& B)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = min(A(i),B(i));
	return tmp;
}

Array_AD max(const Array_AD& A, const Array_AD& B)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = max(A(i),B(i));
	return tmp;
}

Array_AD min(const Array_AD& A, adouble& b)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = min(A(i),b);
	return tmp;
}

Array_AD min(adouble& b, const Array_AD& A)
{
	return min(A,b);
}

Array_AD max(const Array_AD& A, adouble& b)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = max(A(i),b);
	return tmp;
}

Array_AD max(adouble& b, const Array_AD& A)
{
	return max(A,b);
}

adouble sum(const Array_AD& A)
{
	adouble summe=0;
	int    n = A.getSize();
	for (int i=0; i<n; i++)
		summe += A(i);
	return summe;
}


Array_AD abs(const Array_AD& A)
{
	int n = A.getSize();
	Array_AD tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = fabs(A(i));
	return tmp;
}
