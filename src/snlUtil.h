// *** Small Utility Classes and Functions ***

#ifndef SNLUTIL_H
#define SNLUTIL_H

#include "snlPoint.h"
#include "snlVector.h"
#include "snlVersion.h"

#include <iostream>
#include <cmath>

using namespace std;

#ifndef M_PI

	#define M_PI 3.1415926535897932384626433832795

#endif

#ifdef WIN32

	#define isnan _isnan

#endif

const int MAX_BINOMIAL = 64;  // Maximum binomial array dimension.

class binCoefs
{
	// Generate a static array of Binomial Coefficients
	// Array structure is [k][i] where k! / i! ( k - i )!.

	public:

		static int  binCoefArray [ MAX_BINOMIAL ] [ MAX_BINOMIAL ];

		binCoefs();
};

double distToLine ( snlPoint lineStart, snlPoint lineEnd, snlPoint compare );
snlVector projectToLine ( snlPoint lineStart, snlPoint lineEnd, snlPoint compare );

bool isInteriorToTriangle ( snlPoint& testPt, snlPoint& verticeA, snlVector& boundA1, snlVector& boundA2,
							snlPoint& verticeB, snlVector& boundB1, snlVector& boundB2 );

void snlVersion ( int* major, int* minor, int* release );

#endif


