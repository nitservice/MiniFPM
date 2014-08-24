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
 *   Array class adapted for Adol-C
 *
 *************************************************************************/



#ifndef VECTORCLASS_AD_H
#define VECTORCLASS_AD_H

#include <vector>


using namespace std;

class Array_AD
{
private:
	vector<adouble> data;
	int size;

public:
	Array_AD();
	Array_AD(int n);
	Array_AD(adouble& v);
	Array_AD(const Array_AD& V);
	Array_AD(const Array& V);
	Array_AD(int n, adouble* data);
	~Array_AD();


	adouble& operator()(int pos);
	adouble operator()(int pos) const;
	void operator=(const Array_AD& V);
	void operator=(const Array& V);
	void operator=(const adouble& val);
	Array_AD operator+(const Array_AD& V) const;
	Array_AD operator+(const adouble& d) const;
	Array_AD operator-(const Array_AD& V) const;
	Array_AD operator-(const adouble& d) const;
	Array_AD operator*(const Array_AD& V) const;
	Array_AD operator*(const adouble& d) const;
	Array_AD operator/(const Array_AD& V) const;
	Array_AD operator/(const adouble& d) const;
	Array_AD operator<(const Array_AD& V) const;
	Array_AD operator<(adouble& d) const;
	Array_AD operator>(const Array_AD& V) const;
	Array_AD operator>(adouble& d) const;

	void operator<<=(const Array& V);
	void operator>>=(Array& V);


	void remove(int pn);
	void add(int pn, adouble& val);
	void add(int pn);
	void append(adouble& val);
	void setSize(int n);
	void print();


	// Inlines
	int getSize() const {
		return data.size();
	} ;
	adouble value(int pos) const {
		return data[pos];
	};
	bool empty() const {
		return data.empty();
	};
	void clear() {
		data.clear();
	};
	vector<adouble>& vec() {
		return data;
	};

	// Friends
	friend Array_AD operator*(const adouble& d, const Array_AD& V);
	friend Array_AD operator+(const adouble& d, const Array_AD& V);
	friend Array_AD operator-(const adouble& d, const Array_AD& V);
	friend Array_AD operator-(const Array_AD& V);
	friend Array_AD operator/(const adouble& d, const Array_AD& V);

	// Mathematical functions
	friend Array_AD exp(const Array_AD& A);
	friend Array_AD sqrt(const Array_AD& A);
	friend Array_AD sin(const Array_AD& A);
	friend Array_AD cos(const Array_AD& A);
	friend Array_AD pow(const Array_AD& A, int p);
	friend Array_AD abs(const Array_AD& A);
	friend Array_AD min(const Array_AD& A, const Array_AD& B);
	friend Array_AD max(const Array_AD& A, const Array_AD& B);
	friend Array_AD min(const Array_AD& A, adouble& b);
	friend Array_AD max(const Array_AD& A, adouble& b);
	friend Array_AD min(adouble& b, const Array_AD& A);
	friend Array_AD max(adouble& b, const Array_AD& A);
	friend adouble sum(const Array_AD& A);


	friend adouble minval(const Array_AD& vec);
	friend adouble maxval(const Array_AD& vec);
};


#endif
