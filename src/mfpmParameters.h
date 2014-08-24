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
 *   Read and stores parameter given by config.ini
 *
 *************************************************************************/


#ifndef MFPMPARAMETERS_H
#define MFPMPARAMETERS_H

class mfpmParameters
{

public:

	// Parameters
	double AddInner_DeathZone;
	double AddInner_Radius;
	double AddInner_RadiusBnd;
	double AddInner_StartBnd;

	double RemInner_DeathZone;
	double RemInner_Distance;

	double FreeSurf_DeathZone;
	double FreeSurf_RemDist;
	double FreeSurf_AddDist;
	double FreeSurf_InnerDist;
	int FreeSurf_InnerCnt;

	double WallBnd_ConeLength;
	double WallBnd_ConeAngle;

	double InitFixed_Dist;

	double SmoothingLength;
	double SmoothingLengthMin;

	// Constructor
	mfpmParameters();

	// Read config.ini
	void readConfig();

	// Print configuration
	void print();

};


#endif
