// *** NURBS routines that don't belong in classes ***

#include "snlPoint.h"
#include "snlKnotVector.h"
#include "snlSurface.h"
#include "ptrList.h"

#include <iostream>
#include <cmath>

using namespace std;

const int MAX_PROJ_ITER = 64;  // Maximum number of newton iterations for projection functions.

typedef struct
{
	/* surface location */

	knot paramT; // Parameter in t direction.
	knot paramU; // Parameter in u direction.
	snlPoint pt; // Point corresponding to parameters.
	int flag; // Indicates whatever you want ;-)

} sLocn;

typedef struct
{
	/* surface location */

	knot paramT;  // Parameter in t direction.
	knot paramU;  // Parameter in u direction.
	snlPoint pt;  // Point corresponding to parameters.
	int flag;
	basis dist;  // Distance from point to surface.

} projLocn;

bool newtonIterStepSurf(snlPoint* derivPts, snlPoint* projPt, knot* deltaU, knot* deltaV);

bool newtonIterStepCurve(snlPoint* derivPts, snlPoint* projPt, knot* paramDelta);

bool lineIterStepCurve(snlPoint* derivPts, snlPoint* projPt, knot* paramDelta);

int solve2X2LinEqn(double a1, double a2, double a3, double a4, double b1, double b2, double* x1, double* x2);

knot paramDistSqrd(knot t1, knot u1, knot t2, knot u2);

void resolveAmbig(sLocn* projns, int projSize, ptrList <sLocn>* ambig);

ptrList < sLocn >* projPtSurf(snlSurface& surface, snlPoint* toProj, int toProjNum, sLocn* best, double iterTol,
	double cosTol, unsigned maxPass);


