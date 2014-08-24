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
 *   Simple point class storing x/y coordinates and provides simple calculus.
 *
 *************************************************************************/


#ifndef POINTS_H
#define POINTS_H



using namespace std;

class point
{
private:
	void copyData(const point& P);

public:
	double X;
	double Y;
	double Nx;
	double Ny;
	int Bnd;
	int id;
	bool active;
	int meshID;
	double SL;


	// Constructor
	point(int maxNb);
	point(const point& P);
	point();
	~point();

	// Functions
	void calcNormal(const point& P1, const point& P2);
	bool isFreeBnd();
	double dist(const point& P) const;
	double dist2(const point& P) const;

	// Operators
	void operator=(const point& P);
	point operator+(const point& P) const;
	point operator-(const point& P) const;
	point operator*(double lambda) const;
	double operator*(const point& P) const;

	friend point operator*(double lambda, const point& P);

};



#endif

