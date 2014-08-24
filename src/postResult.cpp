#include "postResult.h"
#include <fstream>
#include <vector>

#include "pointSet.h"
#include "scalarField.h"
#include "vectorField.h"


int Suffix_Counter = 0;
int Suffix_Bnd_Counter = 0;
int Suffix_Geo_Counter = 0;

postResult::postResult(pointSet& mesh)
{
	this->mesh = &mesh;
	nbScalarData = 0;
	nbVectorData = 0;
	Suffix_Counter = 0;
	Suffix_Geo_Counter = 0;
	Suffix_Bnd_Counter = 0;
	fileName = "";
	if (sizeof(double)!=8 | sizeof(int)!=4)
	{
		throw "Error postResults: Wrong data format!";
	}
}

postResult::postResult(pointSet& mesh, string fn)
{
	this->mesh = &mesh;
	nbScalarData = 0;
	nbVectorData = 0;
	Suffix_Counter = 0;
	Suffix_Geo_Counter = 0;
	Suffix_Bnd_Counter = 0;
	fileName = fn;
	if (sizeof(double)!=8 | sizeof(int)!=4)
	{
		throw "Error postResults: Wrong data format!";
	}
}

postResult::postResult(pointSet& mesh, int suffix)
{
	this->mesh = &mesh;
	nbScalarData = 0;
	nbVectorData = 0;
	fileName = "";
	Suffix_Counter = suffix;
	Suffix_Geo_Counter = suffix;
	Suffix_Bnd_Counter = suffix;
	if (sizeof(double)!=8 | sizeof(int)!=4)
	{
		throw "Error postResults: Wrong data format!";
	}
}

void postResult::add(scalarField& val)
{
	if (val.getMesh() != mesh)
		throw mfpmExcept(25);
	if (nbScalarData<=nbDataSets)
	{
		sData[nbScalarData] = &val;
	}
	else
	{
		cout << "Error postResults: Too many records!" << endl;
	}
	nbScalarData++;

}

void postResult::add(vectorField& val)
{
	if (val.Xcomp().getMesh() != mesh)
		throw mfpmExcept(25);
	if (nbVectorData<=nbDataSets)
	{
		vData[nbVectorData] = &val;
	}
	else
	{
		cout << "Error postResults: Too many records!" << endl;
	}
	nbVectorData++;

}


void postResult::writeVTK(const std::string& pfileName, bool printMessage)
{
	char ctmp[100];
	double dnull = 0.0;
	double xtmp;


	if (printMessage)
		cout << "Writing file: " << pfileName <<endl;
	ofstream vtk(pfileName.c_str(), ios::binary);
	// Writing header
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "output" << endl;
	vtk << "BINARY" << endl;
	vtk << "DATASET POLYDATA" << endl;
	vtk << "POINTS " << mesh->nbPoints() << " double" << endl;
	// Write points
	unsigned char uctmp[8];
	for (int i=0; i<mesh->nbPoints(); i++)
	{
		xtmp = mesh->x(i);
		convertDoubleChar(uctmp, xtmp);
		vtk.write(reinterpret_cast<char *>(uctmp), 8);
		xtmp = mesh->y(i);
		convertDoubleChar(uctmp, xtmp);
		vtk.write(reinterpret_cast<char *>(uctmp), 8);
		vtk.write(reinterpret_cast<char *>(&dnull), 8);
	}
	vtk << "VERTICES " << mesh->nbPoints() << " " << 2*mesh->nbPoints() << endl;
	unsigned char uitmp[4];
	for (int i=0; i<mesh->nbPoints(); i++)
	{
		convertIntegerChar(uitmp, 1);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, i);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
	}
	// Write scalar data
	if (nbScalarData>0 | nbVectorData>0)
		vtk << "POINT_DATA " << mesh->nbPoints() << endl;
	for (int scalDat=0; scalDat<nbScalarData; scalDat++)
	{
		vtk << "SCALARS " << sData[scalDat]->name << " double 1" << endl;
		vtk << "LOOKUP_TABLE default" << endl;
		for (int i=0; i<mesh->nbPoints(); i++)
		{
			convertDoubleChar(uctmp, sData[scalDat]->value(i));
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
		}
	}
	// Write Array data
	for (int vecDat=0; vecDat<nbVectorData; vecDat++)
	{
		vtk << "VECTORS " << vData[vecDat]->name << " double" << endl;
		scalarField& valX = vData[vecDat]->Xcomp();
		scalarField& valY = vData[vecDat]->Ycomp();
		for (int i=0; i<mesh->nbPoints(); i++)
		{
			convertDoubleChar(uctmp, valX(i));
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
			convertDoubleChar(uctmp, valY(i));
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
			convertDoubleChar(uctmp, 0.0);
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
		}
	}
	vtk.close();
}

void postResult::writeVTK(const char *cfileName, bool printMessage)
{
	writeVTK(string(cfileName), printMessage);
}


void postResult::writeVTK(int suffix, bool printMessage)
{
	char buf[100];
	sprintf(buf, "%s_%5.5i.vtk", fileName.c_str(), suffix );
	writeVTK(buf, printMessage);
}


void postResult::writeVTK(bool printMessage)
{
	char buf[100];
	sprintf(buf, "%s_%5.5i.vtk", fileName.c_str(), Suffix_Counter );
	writeVTK(buf, printMessage);
	Suffix_Counter++;
}


void postResult::writeGeometryVTK(const string& pfileName, bool printMsg)
{
	char ctmp[100];
	double dnull = 0.0;
	double xtmp;

	if (printMsg)
		cout << "Writing file: " << pfileName <<endl;
	ofstream vtk(pfileName.c_str(), ios::binary);
	// Writing header
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "output" << endl;
	vtk << "BINARY" << endl;
	vtk << "DATASET POLYDATA" << endl;
	vtk << "POINTS " << mesh->nbPoints() << " double" << endl;
	// Write points
	unsigned char uctmp[8];
	for (int i=0; i<mesh->nbPoints(); i++)
	{
		xtmp = mesh->x(i);
		convertDoubleChar(uctmp, xtmp);
		vtk.write(reinterpret_cast<char *>(uctmp), 8);
		xtmp = mesh->y(i);
		convertDoubleChar(uctmp, xtmp);
		vtk.write(reinterpret_cast<char *>(uctmp), 8);
		vtk.write(reinterpret_cast<char *>(&dnull), 8);
	}
	vtk << "VERTICES " << mesh->nbPoints() << " " << 2*mesh->nbPoints() << endl;
	unsigned char uitmp[4];
	for (int i=0; i<mesh->nbPoints(); i++)
	{
		convertIntegerChar(uitmp, 1);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, i);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
	}
	vtk << "LINES " << 2*mesh->getNbBndPoints() << " " << 2*3*mesh->getNbBndPoints() << endl;
	int rnb,lnb;
	for (int i=0; i<mesh->getNbBndPoints(); i++)
	{
		convertIntegerChar(uitmp, 2);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		mesh->findLeftRightBndPoint(i, rnb, lnb);
		convertIntegerChar(uitmp, i);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, rnb);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);

		convertIntegerChar(uitmp, 2);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, i);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, lnb);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
	}
	// Write scalar data
	vtk << "POINT_DATA " << mesh->nbPoints() << endl;
	vtk << "SCALARS BndID double 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl;
	for (int i=0; i<mesh->nbPoints(); i++)
	{
		convertDoubleChar(uctmp, mesh->getBND(i));
		vtk.write(reinterpret_cast<char *>(&uctmp), 8);
	}

	// Write Array data
	vtk << "VECTORS Normals double" << endl;
	for (int i=0; i<mesh->nbPoints(); i++)
	{
		convertDoubleChar(uctmp, mesh->getNormal(0,i));
		vtk.write(reinterpret_cast<char *>(&uctmp), 8);
		convertDoubleChar(uctmp, mesh->getNormal(1,i));
		vtk.write(reinterpret_cast<char *>(&uctmp), 8);
		convertDoubleChar(uctmp, 0.0);
		vtk.write(reinterpret_cast<char *>(&uctmp), 8);
	}

	vtk.close();
}

void postResult::writeGeometryVTK(const char *cfileName, bool printMessage)
{
	writeGeometryVTK(string(cfileName),printMessage);
}


void postResult::writeGeometryVTK(int suffix, bool printMessage)
{
	char buf[100];
	sprintf(buf, "%s_geo_%5.5i.vtk", fileName.c_str(), suffix );
	writeVTK(buf, printMessage);
}


void postResult::writeGeometryVTK(bool printMessage)
{
	char buf[100];
	sprintf(buf, "%s_geo_%5.5i.vtk", fileName.c_str(), Suffix_Geo_Counter );
	writeVTK(buf, printMessage);
	Suffix_Geo_Counter++;
}



void postResult::writeBoundaryVTK(const string& pfileName, bool printMsg)
{
	char ctmp[100];
	double dnull = 0.0;
	double xtmp;

	if (printMsg)
		cout << "Writing file: " << pfileName <<endl;
	ofstream vtk(pfileName.c_str(), ios::binary);
	// Writing header
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "output" << endl;
	vtk << "BINARY" << endl;
	vtk << "DATASET POLYDATA" << endl;
	vtk << "POINTS " << mesh->getNbBndPoints() << " double" << endl;
	// Write points
	unsigned char uctmp[8];
	for (int i=0; i<mesh->getNbBndPoints(); i++)
	{
		xtmp = mesh->x(i);
		convertDoubleChar(uctmp, xtmp);
		vtk.write(reinterpret_cast<char *>(uctmp), 8);
		xtmp = mesh->y(i);
		convertDoubleChar(uctmp, xtmp);
		vtk.write(reinterpret_cast<char *>(uctmp), 8);
		vtk.write(reinterpret_cast<char *>(&dnull), 8);
	}
	vtk << "LINES " << 2*mesh->getNbBndPoints() << " " << 2*3*mesh->getNbBndPoints() << endl;
	unsigned char uitmp[4];
	int rnb,lnb;
	for (int i=0; i<mesh->getNbBndPoints(); i++)
	{
		convertIntegerChar(uitmp, 2);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		mesh->findLeftRightBndPoint(i, rnb, lnb);
		convertIntegerChar(uitmp, i);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, rnb);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);

		convertIntegerChar(uitmp, 2);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, i);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
		convertIntegerChar(uitmp, lnb);
		vtk.write(reinterpret_cast<char *>(uitmp), 4);
	}
	// Write scalar data
	if (nbScalarData>0 | nbVectorData>0)
		vtk << "POINT_DATA " << mesh->getNbBndPoints() << endl;
	for (int scalDat=0; scalDat<nbScalarData; scalDat++)
	{
		vtk << "SCALARS " << sData[scalDat]->name << " double 1" << endl;
		vtk << "LOOKUP_TABLE default" << endl;
		for (int i=0; i<mesh->getNbBndPoints(); i++)
		{
			convertDoubleChar(uctmp, sData[scalDat]->value(i));
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
		}
	}
	// Write Array data
	for (int vecDat=0; vecDat<nbVectorData; vecDat++)
	{
		vtk << "VECTORS " << vData[vecDat]->name << " double" << endl;
		scalarField& valX = vData[vecDat]->Xcomp();
		scalarField& valY = vData[vecDat]->Ycomp();
		for (int i=0; i<mesh->getNbBndPoints(); i++)
		{
			convertDoubleChar(uctmp, valX(i));
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
			convertDoubleChar(uctmp, valY(i));
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
			convertDoubleChar(uctmp, 0.0);
			vtk.write(reinterpret_cast<char *>(&uctmp), 8);
		}
	}

	vtk.close();
}

void postResult::writeBoundaryVTK(const char *cfileName, bool printMessage)
{
	writeBoundaryVTK(string(cfileName), printMessage);
}



void postResult::writeBoundaryVTK(int suffix, bool printMessage)
{
	char buf[100];
	sprintf(buf, "%s_bnd_%5.5i.vtk", fileName.c_str(), suffix );
	writeBoundaryVTK(buf, printMessage);
}


void postResult::writeBoundaryVTK(bool printMessage)
{
	char buf[100];
	sprintf(buf, "%s_bnd_%5.5i.vtk", fileName.c_str(), Suffix_Bnd_Counter );
	writeBoundaryVTK(buf, printMessage);
	Suffix_Bnd_Counter++;
}



void postResult::convertDoubleChar(unsigned char co[8], double di)
{
	unsigned char* xtm = (unsigned char*)&di;
	for (int i=0; i<8; i++)
	{
		co[i] = xtm[7-i];
	}
}


void postResult::convertIntegerChar(unsigned char co[4], int di)
{
	unsigned char* xtm = (unsigned char*)&di;
	for (int i=0; i<4; i++)
	{
		co[i] = xtm[3-i];
	}

}


void postResult::setSuffix(int s1, int s2, int s3)
{
	if ( s1 >= 0) Suffix_Counter = s1;
	if ( s2 >= 0) Suffix_Bnd_Counter = s2;
	if ( s3 >= 0) Suffix_Geo_Counter = s3;
}


void postResult::setSuffix(int s1)
{
	Suffix_Counter = s1;
	Suffix_Bnd_Counter = s1;
	Suffix_Geo_Counter = s1;
}

