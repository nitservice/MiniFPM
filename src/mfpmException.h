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
 *   Exception class for MiniFPM-Exceptions
 *
 *************************************************************************/


#ifndef MFPMEXCEPTIONS_H
#define MFPMEXCEPTIONS_H

#define MFPM_EXCEPT __FILE__, __LINE__, __FUNCTION__
#define mfpmExcept(N) mfpmException(MFPM_EXCEPT, N)


#include <iostream>
#include <exception>

using namespace std;

class mfpmException: public exception
{
private:
	string errorMsg;
	int errCode;

public:
	mfpmException(string fn, int ln, string funcN, int exNb);
	virtual ~mfpmException() throw();

	virtual const char* what() const throw();

	const int error() const throw();

};

#endif
