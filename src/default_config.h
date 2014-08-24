/*****************************************************************
* Configuration file for pointSet.cpp                            *
* ===================================                            *
*                                                                *
* This is the default setting if no config.ini is present.       *
*                                                                *
* All values are: ... * SmoothingLength                          *
*                                                                *
******************************************************************/

// Adapt inner points
#define __AddInner_DeathZone   0.15	 	// Defines tolerance for point in/out (inner point <-> wall) [0.15]
#define __AddInner_Radius      0.4		// Radius of ball for hole-check [0.4]
#define __AddInner_RadiusBnd   0.2		// Radius of ball near boundary [???]
#define __AddInner_StartBnd    2.0		// Start range (...*SL) for __RadiusBnd (linear decrease)

#define __RemInner_DeathZone   0.15 	// Defines tolerance for point in/out (inner point <-> wall) [0.15]
#define __RemInner_Distance    0.2		// Points closer than ... are removed [0.2]


// Adapt free surface points
#define __FreeSurf_DeathZone   0.2		// Defines tolerance for point in/out (free surface) [0.1]
#define __FreeSurf_RemDist     0.2		// Distance where free surface points are removed [0.15]
#define __FreeSurf_AddDist     0.35		// Distance where free surface points are added (should be >__FreeSurf_RemDist) [0.3]
#define __FreeSurf_InnerDist   0.8		// Max dist. to inner point
#define __FreeSurf_InnerCnt    2		// Inner points needed in range of __InnerDist


// Acitvate / Deactivate geometric points
#define __WallBnd_ConeLength   1.0
#define __WallBnd_ConeAngle    5.0


// Initialisation
#define __InitFixed_Dist       0.35		// Distance of "wall" boundary points [0.3 / 0.4]


