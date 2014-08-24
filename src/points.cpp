#include "points.h"
#include <iostream>
#include <vector>
#include <cmath>


point::point(int maxNb)
{
	X = 0;
	Y = 0;
	Nx = 0;
	Ny = 0;
	Bnd = 0;
	id = 0;
	SL = 0;
	active = true;
}


point::point()
{
	X = 0;
	Y = 0;
	Nx = 0;
	Ny = 0;
	Bnd = 0;
	id = 0;
	meshID = -1;
	SL = 0;
	active = true;
}


point::point(const point& P)
{
	copyData(P);
}


point::~point()
{

}



void point::copyData(const point& P)
{
	X = P.X;
	Y = P.Y;
	Nx = P.Nx;
	Ny = P.Ny;
	Bnd = P.Bnd;
	id = P.id;
	SL = P.SL;
	active = P.active;

}


void point::operator=(const point& P)
{
	copyData(P);
}


void point::calcNormal(const point& P1, const point& P2)
{
	// Central
	Nx = P2.Y - P1.Y;
	Ny = P1.X - P2.X;

	double length = sqrt(Nx*Nx + Ny*Ny);
	Nx = Nx / length;
	Ny = Ny / length;

}




double point::dist(const point& P) const
{
	return sqrt( pow(P.X-X,2) + pow(P.Y-Y,2) );
}


double point::dist2(const point& P) const
{
	return pow(P.X-X,2) + pow(P.Y-Y,2);
}



bool point::isFreeBnd()
{
	if (Bnd<200)
		return false;
	else
		return true;
}



point point::operator+(const point& P) const
{
	point pTmp;
	pTmp.X = X + P.X;
	pTmp.Y = Y + P.Y;
	return pTmp;
}

point point::operator-(const point& P) const
{
	point pTmp;
	pTmp.X = X - P.X;
	pTmp.Y = Y - P.Y;
	return pTmp;
}

point point::operator*(double lambda) const
{
	point pTmp;
	pTmp.X = lambda*X;
	pTmp.Y = lambda*Y;
	return pTmp;
}

double point::operator*(const point& P) const
{
	return P.X*X + P.Y*Y;
}


point operator*(double lambda, const point& P)
{
	return P*lambda;
}



