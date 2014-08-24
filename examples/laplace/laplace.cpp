#include <cmath>
#include <iostream>
#include <cstdio>

#include "pointSet.h"
#include "scalarField.h"
#include "postResult.h"
#include "primitive_operators.h"
#include "Array.h"
#include "mfpmMatrix.h"
#include "operators.h"

#include "matrixUmfpack.h"
#include "mfpmParameters.h"

using namespace std;

int main(int argc, char** argv)
{
	
	pointSet mesh("geometry.msh");	
	postResult post(mesh);

	scalarField u(mesh, "u");
	post.add(u);

	mfpmMatrix A = -1.0*laplace(mesh) == 1.0;
	A.setDirichlet(10, 0.0);
	A.setDirichlet(20, 0.0);
	A.setDirichlet(30, 0.0);
	A.setDirichlet(40, 0.0);
	A.solveUMFPACK(u);

	post.writeGeometryVTK("result_geo.vtk");
	post.writeVTK("result_data.vtk");


	return 0;

}



