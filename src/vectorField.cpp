#include "vectorField.h"
#include "scalarField.h"
#include "Array.h"

vectorField::vectorField()
{
	Xval = new scalarField();
	Yval = new scalarField();
	name = "tmpVec";
}

vectorField::vectorField(pointSet& mesh)
{
	Xval = new scalarField(mesh, "tmpVecX");
	Yval = new scalarField(mesh, "tmpVecY");
	name = "tmpVec";
}

vectorField::vectorField(pointSet& mesh, const string& name)
{
	Xval = new scalarField(mesh, name+"_X");
	Yval = new scalarField(mesh, name+"_Y");
	this->name = name;
}


vectorField::~vectorField()
{
	delete Xval;
	delete Yval;
}


double& vectorField::Xcomp(int i)
{
	return Xval->operator()(i);
}

double& vectorField::Ycomp(int i)
{
	return Yval->operator()(i);
}

