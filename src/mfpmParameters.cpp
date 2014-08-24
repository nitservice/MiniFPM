#include "mfpmParameters.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;



mfpmParameters::mfpmParameters()
{
	AddInner_DeathZone = 0.15;
	AddInner_Radius    = 0.40;
	AddInner_RadiusBnd = 0.20;
	AddInner_StartBnd  = 2.0;

	RemInner_DeathZone = 0.15;
	RemInner_Distance  = 0.2;

	FreeSurf_DeathZone = 0.2;
	FreeSurf_RemDist   = 0.2;
	FreeSurf_AddDist   = 0.35;
	FreeSurf_InnerDist = 0.8;
	FreeSurf_InnerCnt  = 2;

	WallBnd_ConeLength = 1.0;
	WallBnd_ConeAngle  = 5.0;

	InitFixed_Dist     = 0.35;

	SmoothingLength    = 0.2;

	SmoothingLengthMin = -1;

	readConfig();

}


void mfpmParameters::readConfig()
{
	ifstream conf("config.ini");
	if (conf.fail())
	{
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "!! Configuration file not found! Using default values!" << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		return;
	}

	char cbuf[100];
	string sbuf;
	int pos, lcnt = 0;
	double inVal;
	while (!conf.eof())
	{
		lcnt++;
		conf.getline(cbuf,100);
		int n = conf.gcount();
		if (n<=2)
			continue;
		if (cbuf[0] == '#')
			continue;
		sbuf = cbuf;
		pos = sbuf.find("=");
		if (pos < 0)
		{
			cout << "mfpmParameters::readConfig: invalid entry found at line " << lcnt << endl;
			continue;
		}
		inVal = atof(sbuf.substr(pos+1).c_str());
		if (sbuf.find("AddInner_DeathZone")!=string::npos)
			AddInner_DeathZone = inVal;
		else if (sbuf.find("AddInner_RadiusBnd")!=string::npos)
			AddInner_RadiusBnd = inVal;
		else if (sbuf.find("AddInner_Radius")!=string::npos)
			AddInner_Radius = inVal;
		else if (sbuf.find("AddInner_StartBnd")!=string::npos)
			AddInner_StartBnd = inVal;
		else if (sbuf.find("RemInner_DeathZone")!=string::npos)
			RemInner_DeathZone = inVal;
		else if (sbuf.find("RemInner_Distance")!=string::npos)
			RemInner_Distance = inVal;
		else if (sbuf.find("FreeSurf_DeathZone")!=string::npos)
			FreeSurf_DeathZone = inVal;
		else if (sbuf.find("FreeSurf_RemDist")!=string::npos)
			FreeSurf_RemDist = inVal;
		else if (sbuf.find("FreeSurf_AddDist")!=string::npos)
			FreeSurf_AddDist = inVal;
		else if (sbuf.find("FreeSurf_InnerDist")!=string::npos)
			FreeSurf_InnerDist = inVal;
		else if (sbuf.find("FreeSurf_InnerCnt")!=string::npos)
			FreeSurf_InnerCnt = (int)inVal;
		else if (sbuf.find("WallBnd_ConeLength")!=string::npos)
			WallBnd_ConeLength = inVal;
		else if (sbuf.find("WallBnd_ConeAngle")!=string::npos)
			WallBnd_ConeAngle = inVal;
		else if (sbuf.find("InitFixed_Dist")!=string::npos)
			InitFixed_Dist = inVal;
		else if (sbuf.find("SmoothingLength")!=string::npos)
			SmoothingLength = inVal;
		else if (sbuf.find("hmax")!=string::npos)
			SmoothingLength = inVal;
		else if (sbuf.find("hmin")!=string::npos)
			SmoothingLengthMin = inVal;


	}

	if (SmoothingLengthMin < 0)
		SmoothingLengthMin = SmoothingLength;


}



void mfpmParameters::print()
{
	cout << "Parameters used:" << endl;
	cout << "================" << endl;
	cout << "SmoothingLength:       " << SmoothingLength << endl;
	cout << "AddInner_DeathZone:    " << AddInner_DeathZone << endl;
	cout << "AddInner_Radius:       " << AddInner_Radius << endl;
	cout << "AddInner_RadiusBnd:    " << AddInner_RadiusBnd << endl;
	cout << "AddInner_StartBnd:     " << AddInner_StartBnd << endl;
	cout << "RemInner_DeathZone:    " << RemInner_DeathZone << endl;
	cout << "RemInner_Distance:     " << RemInner_Distance << endl;
	cout << "FreeSurf_DeathZone:    " << FreeSurf_DeathZone << endl;
	cout << "FreeSurf_RemDist:      " << FreeSurf_RemDist << endl;
	cout << "FreeSurf_AddDist:      " << FreeSurf_AddDist << endl;
	cout << "FreeSurf_InnerDist:    " << FreeSurf_InnerDist << endl;
	cout << "FreeSurf_InnerCnt:     " << FreeSurf_InnerCnt << endl;
	cout << "WallBnd_ConeLength:    " << WallBnd_ConeLength << endl;
	cout << "WallBnd_ConeAngle:     " << WallBnd_ConeAngle << endl;
	cout << "InitFixed_Dist:        " << InitFixed_Dist << endl;

}
