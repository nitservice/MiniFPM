#include "mfpmException.h"
#include <cstdio>


mfpmException::mfpmException(string fn, int ln, string funcN, int exNb)
{
	string msg;
	errCode = exNb;
	switch (exNb)
	{
	case 0:
		msg = "No boundary point";
		break;
	case 1:
		msg = "No inner point";
		break;
	case 2:
		msg = "Index out of range";
		break;
	case 3:
		msg = "Invalid component index";
		break;
	case 4:
		msg = "Requested point not found in neighbourhood";
		break;
	case 5:
		msg = "Geometry point not active";
		break;
	case 6:
		msg = "Point ouside of bounding box";
		break;
	case 10:
		msg = "Field not found";
		break;
	case 20:
		msg = "Not enough points present";
		break;
	case 21:
		msg = "Not enough boundary points present";
		break;
	case 22:
		msg = "Too many points added";
		break;
	case 23:
		msg = "Smoothing length too high";
		break;
	case 24:
		msg = "Cannot assign same mesh. 'M = M' not possible";
		break;
	case 25:
		msg = "PointSets differ";
		break;
	case 26:
		msg = "Order of interpolation not high enough";
		break;
	case 27:
		msg = "Unkown error solving least squares problem";
		break;
	case 100:
		msg = "Wrong GMSH format of file not found";
		break;
	case 200:
		msg = "Umfpack matrix sizes differ";
		break;

	default:
		msg = "Unknown error";
	}

	char buf[200];
	sprintf(buf, "Error %s(%i):: %s: %s!", fn.c_str(), ln, funcN.c_str(), msg.c_str());
	errorMsg = buf;
}

mfpmException::~mfpmException() throw()
{
}

const char* mfpmException::what() const throw()
{
	return errorMsg.c_str();
}

const int mfpmException::error() const throw()
{
	return errCode;
}


