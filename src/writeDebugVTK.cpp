#include "writeDebugVTK.h"

#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

int counter = 0;

void writeDebugVTK(int np, double* x, double* y, double* val)
{
	char buf[100];
	sprintf(buf, "debug_%4.4i.vtk", counter++);
	ofstream vtk(buf);
	// Writing header
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "output" << endl;
	vtk << "ASCII" << endl;
	vtk << "DATASET POLYDATA" << endl;
	vtk << "POINTS " << np << " double" << endl;

	// Point data
	for (int i=0; i<np; i++)
	{
		sprintf(buf, "%12.4E %12.4E %12.4E", x[i], y[i], 0.0);
		vtk << buf << endl;
		//vtk << x[i] << " " << y[i] << " " << 0.0 << endl;
	}
	// Vertices
	vtk << "VERTICES " << np << " " << 2*np << endl;
	for (int i=0; i<np; i++)
	{
		vtk << "1 " << i+1 << endl;
	}

	// Point data
	vtk << "POINT_DATA " << np << endl;
	vtk << "SCALARS debug double 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl;
	for (int i=0; i<np; i++)
	{
		sprintf(buf, "%12.4E", val[i]);
		vtk << buf << endl;
	}
}



void writeDebugVTK(int np, double* x, double* y)
{
	char buf[100];
	sprintf(buf, "debug_%4.4i.vtk", counter++);
	ofstream vtk(buf);
	// Writing header
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "output" << endl;
	vtk << "ASCII" << endl;
	vtk << "DATASET POLYDATA" << endl;
	vtk << "POINTS " << np << " double" << endl;

	// Point data
	for (int i=0; i<np; i++)
	{
		sprintf(buf, "%12.4E %12.4E %12.4E", x[i], y[i], 0.0);
		vtk << buf << endl;
		//vtk << x[i] << " " << y[i] << " " << 0.0 << endl;
	}
	// Vertices
	vtk << "VERTICES " << np << " " << 2*np << endl;
	for (int i=0; i<np; i++)
	{
		vtk << "1 " << i+1 << endl;
	}

}




void writeDebugVTKwithLine(int np, double* x, double* y)
{
	char buf[100];
	sprintf(buf, "debug_%4.4i.vtk", counter++);
	ofstream vtk(buf);
	// Writing header
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "output" << endl;
	vtk << "ASCII" << endl;
	vtk << "DATASET POLYDATA" << endl;
	vtk << "POINTS " << np << " double" << endl;

	// Point data
	for (int i=0; i<np; i++)
	{
		sprintf(buf, "%12.4E %12.4E %12.4E", x[i], y[i], 0.0);
		vtk << buf << endl;
		//vtk << x[i] << " " << y[i] << " " << 0.0 << endl;
	}
	// Vertices
	vtk << "LINES " << np << " " << 3*np << endl;
	for (int i=0; i<np-1; i++)
	{
		vtk << "2 " << i << " " << i+1 << endl;
	}
	vtk << "2 " << np-1 << " " << 0 << endl;

}

