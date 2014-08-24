#include "Array.h"
#include <cmath>
#include <iostream>

Array::Array(int n)
{
	size = n;
	data.resize(n,0.0);
}

Array::Array(double v)
{
	cout << "Warning: Array constructor with double value! -> empty Array has been created!" << endl;
	size = 0;
}

Array::Array(const Array& V)
{
	size = V.getSize();
	data.reserve(size);
	for ( int i=0; i<size; i++)
		data.push_back(V.value(i));
}

Array::Array(int n, double* val)
{
	size = n;
	data.reserve(size);
	for ( int i=0; i<size; i++)
		data.push_back(val[i]);
}

Array::Array()
{
	size = 0;
}

Array::~Array()
{
	data.clear();
}


double& Array::operator()(int pos)
{
	return data[pos];
}

double Array::operator()(int pos) const
{
	return data[pos];
}



void Array::operator=(const Array& V)
{
	if (data.empty())
	{
		size = V.getSize();
		data.reserve(size);
		for ( int i=0; i<size; i++)
			data.push_back(V.value(i));
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

void Array::operator=(const double val)
{
	if (data.empty())
	{
		throw "Error Array:: operator=: object not initialised!";
	}
	else
	{
		for (int i=0; i<size; i++)
			data[i] = val;
	}
}


Array Array::operator+(const Array& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array class:: operator+: sizes differ";
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] + V.value(i);
	return A;
}


Array Array::operator+(const double v) const
{
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] + v;
	return A;
}

Array operator+(const double d, const Array& V)
{
	return V+d;
}


Array Array::operator-(const Array& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array class:: operator+: sizes differ";
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] - V.value(i);
	return A;
}

Array Array::operator-(const double v) const
{
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] - v;
	return A;
}

Array operator-(const double d, const Array& V)
{
	Array A(V.getSize());
	for (int i=0; i<V.getSize(); i++)
		A(i) = d - V.value(i);
	return A;
}


Array operator-(const Array& V)
{
	Array A(V.getSize());
	for (int i=0; i<V.getSize(); i++)
		A(i) = - V.value(i);
	return A;
}



Array Array::operator*(const Array& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array class:: operator*: sizes differ";
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] * V.value(i);
	return A;
}


Array Array::operator*(const double V) const
{
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] * V;
	return A;
}

Array operator*(const double d, const Array& V)
{
	return V*d;
}



Array Array::operator/(const Array& V) const
{
	if (this->size!=V.getSize())
		throw "Error Array class:: operator/: sizes differ";
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] / V.value(i);
	return A;
}

Array Array::operator/(const double d) const
{
	Array A(this->size);
	for (int i=0; i<this->size; i++)
		A(i) = this->data[i] / d;
	return A;
}


Array operator/(const double d, const Array& V)
{
	Array A(V.getSize());
	for (int i=0; i<V.getSize(); i++)
		A(i) = d / V.value(i);
	return A;
}

Array Array::operator<(const Array& A) const
{
	Array tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] < A(i) )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}


Array Array::operator<(double d) const
{
	Array tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] < d )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}


Array Array::operator>(const Array& A) const
{
	Array tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] > A(i) )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}


Array Array::operator>(double d) const
{
	Array tmp(this->size);
	for (int i=0; i<this->size; i++)
		if ( this->data[i] > d )
			tmp(i) = 1.0;
		else
			tmp(i) = 0.0;
	return tmp;
}





void Array::setSize(int n)
{
	size = n;
	data.resize(n, 0.0);
}


void Array::remove(int pn)
{
	if (pn<size)
	{
		vector<double>::iterator it;
		it = data.begin() + pn;
		data.erase(it);
		size=data.size();
	}
	else
		throw "Error Array::remove: invalid point number";
}



void Array::add(int pn, double val)
{
	if (pn<=size & pn>=0)
	{
		data.insert(data.begin() + pn, val);
		size=data.size();
	}
	else
		throw "Error Array::add(int, double): invalid point number";
}


void Array::add(int pn)
{
	if (pn<=size & pn>=0)
	{
		data.insert(data.begin() + pn, 0.0);
		size=data.size();
	}
	else
		throw "Error Array::add(int, double): invalid point number";
}


void Array::append(double val)
{
	data.push_back(val);
	size=data.size();
}


double minval(const Array& vec)
{
	double mv = vec(0);
	for (int i=1; i<vec.size; i++)
		if (mv>vec(i))
			mv = vec(i);
	return mv;
}

double maxval(const Array& vec)
{
	double mv = vec(0);
	for (int i=1; i<vec.size; i++)
		if (mv<vec(i))
			mv = vec(i);
	return mv;
}

void Array::print()
{
	int n = getSize();
	cout << "Array values: " << n << endl;
	cout << "============================================" << endl;
	for (int i=0; i<n; i++)
		cout << data[i] << endl;
	cout << "============================================" << endl;
}


void Array::vec(double* arr)
{
	for(int i=0; i<size; i++)
		arr[i] = data[i];
}


/********************************
 * Math operators, e.g. exp, pow
 ********************************/
Array exp(const Array& A)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = exp(A(i));
	return tmp;
}

Array sqrt(const Array& A)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = sqrt(A(i));
	return tmp;
}

Array atan(const Array& A)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = atan(A(i));
	return tmp;
}

Array sin(const Array& A)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = sin(A(i));
	return tmp;
}

Array cos(const Array& A)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = cos(A(i));
	return tmp;
}

Array pow(const Array& A, int p)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = pow(A(i), p);
	return tmp;
}



Array min(const Array& A, const Array& B)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = min(A(i),B(i));
	return tmp;
}

Array max(const Array& A, const Array& B)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = max(A(i),B(i));
	return tmp;
}

Array min(const Array& A, double b)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = min(A(i),b);
	return tmp;
}

Array min(double b, const Array& A)
{
	return min(A,b);
}

Array max(const Array& A, double b)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = max(A(i),b);
	return tmp;
}

Array max(double b, const Array& A)
{
	return max(A,b);
}

double sum(const Array& A)
{
	double summe=0;
	int    n = A.getSize();
	for (int i=0; i<n; i++)
		summe += A(i);
	return summe;
}


Array abs(const Array& A)
{
	int n = A.getSize();
	Array tmp(n);
	for (int i=0; i<n; i++)
		tmp(i) = abs(A(i));
	return tmp;
}
