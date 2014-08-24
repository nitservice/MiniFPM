/*************************************************************************
 *
 *  This file is part of MiniFPM
 *  Copyright (C) 2014 Jan Marburger
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Description:
 *   Array class is a container for data supplying mathematical operations.
 *
 *************************************************************************/


#ifndef VECTORCLASS_H
#define VECTORCLASS_H

#include <vector>


using namespace std;

class Array
{
private:
	Array(double n);
	vector<double> data;
	int size;

public:
	Array();
	Array(int n);
	Array(const Array& V);
	Array(int n, double* data);
	~Array();


	double& operator()(int pos);
	double operator()(int pos) const;
	void operator=(const Array& V);
	void operator=(const double val);
	Array operator+(const Array& V) const;
	Array operator+(const double d) const;
	Array operator-(const Array& V) const;
	Array operator-(const double d) const;
	Array operator*(const Array& V) const;
	Array operator*(const double d) const;
	Array operator/(const Array& V) const;
	Array operator/(const double d) const;
	Array operator<(const Array& V) const;
	Array operator<(double d) const;
	Array operator>(const Array& V) const;
	Array operator>(double d) const;


	void remove(int pn);
	void add(int pn, double val);
	void add(int pn);
	void append(double val);
	void setSize(int n);
	void print();


	// Inlines
	int getSize() const {
		return data.size();
	} ;
	double value(int pos) const {
		return data[pos];
	};
	bool empty() const {
		return data.empty();
	};
	void clear() {
		data.clear();
	};
	vector<double>& vec() {
		return data;
	};
	void vec(double* arr);

	// Friends
	friend Array operator*(const double d, const Array& V);
	friend Array operator+(const double d, const Array& V);
	friend Array operator-(const double d, const Array& V);
	friend Array operator-(const Array& V);
	friend Array operator/(const double d, const Array& V);

	// Mathematical functions
	friend Array exp(const Array& A);
	friend Array sqrt(const Array& A);
	friend Array sin(const Array& A);
	friend Array cos(const Array& A);
	friend Array atan(const Array& A);
	friend Array pow(const Array& A, int p);
	friend Array abs(const Array& A);
	friend Array min(const Array& A, const Array& B);
	friend Array max(const Array& A, const Array& B);
	friend Array min(const Array& A, double b);
	friend Array max(const Array& A, double b);
	friend Array min(double b, const Array& A);
	friend Array max(double b, const Array& A);
	friend double sum(const Array& A);


	friend double minval(const Array& vec);
	friend double maxval(const Array& vec);
};


#endif
