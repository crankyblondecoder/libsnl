// *** General NURBS Curve ***

#include "snlCurve.h"
#include "snlUtil.h"
#include "snlSquareLinear.h"

snlCurve::~snlCurve()
{
	if(ctrlPtNet) delete ctrlPtNet;
	if(knotVect) delete knotVect;
}

snlCurve::snlCurve()
{
	deg = 0;

	ctrlPtNet = 0;
	knotVect = 0;
}

snlCurve::snlCurve(const snlCurve& copyFrom)
{
	// Copy constructor.
	// -----------------

	ctrlPtNet = new snlCtrlPointNetCurve(*(copyFrom.ctrlPtNet));

	knotVect = new snlKnotVector(*(copyFrom.knotVect));

	deg = copyFrom.deg;
}

snlCurve::snlCurve(int degree, unsigned size, snlPoint& start, snlPoint& end)
{
	// Generate new curve.
	// -------------------
	// deg:    Degree of curve.
	// size:   Number of control points curve has.
	// start:  Starting point of curve.
	// end:    Ending point of curve.
	//
	// Notes:  Draws a straight line.

	deg = degree;

	ctrlPtNet = new snlCtrlPointNetCurve(size, start, end);

	knotVect = new snlKnotVector(0.0, 1.0, size + degree + 1, degree);
}

snlCurve::snlCurve(int degree, unsigned size, snlCtrlPoint* points, knot* knots)
{
	// Generate new curve.
	// -------------------
	// deg:        Degree of curve.
	// size:       Number of control points curve has.
	// points:     Points to use.
	// knots:      Knots to use.
	//
	// Notes:      Does not copy points and knots. So don't delete them elsewhere.

	deg = degree;

	ctrlPtNet = new snlCtrlPointNetCurve(points, size, false);

	if(knots)
		knotVect = new snlKnotVector(knots, size + degree + 1, degree);
	else
		knotVect = new snlKnotVector(0.0, 1.0, size + degree + 1, degree);
}

snlCurve::snlCurve(unsigned size, snlCtrlPoint* points, snlKnotVector* knotVector)
{
	// Generate new curve.
	// -------------------
	// size:        Size of control point array.
	// points:      Control point array.
	// knotVect:    Knot vector to use.
	//
	// Notes:       Does not copy points and knot vector. So don't delete them elsewhere.

	deg = knotVector -> getDegree();

	knotVect = knotVector;

	ctrlPtNet = new snlCtrlPointNetCurve(points, size, false);
}

snlCurve::snlCurve(snlPoint* points, unsigned size, int fittingType, int degree, bool closedLoop, knot** retParams)
{
	// Interpolated / approximated curve.
	// ----------------------------------
	// points:      Points to interpolate between.
	// size:        Number of points.
	// fittingType: Type of interpolation from SNL_FITTING_TYPES.
	// degree:      Resultant curve should be this degree.
	// closedLoop:  The points specify a closed loop that should join smoothly.
	// retParams:   Pointer to pointer that points to array of parameters that correspond to given points.
	//
	// Notes:       Array returned via retParams should be deleted by calling function.

	ctrlPtNet = 0;
	knotVect = 0;

	switch(fittingType)
	{
		case SNL_GLOBAL_INTERPOLATION:
		case SNL_GLOBAL_INTERP_CENTRIFUGAL:

			if(closedLoop)
				globalInterpClosedLoop(fittingType, points, size, degree, retParams);
			else
				globalInterp(fittingType, points, size, degree, retParams);

			break;

		case SNL_LOCAL_INTERPOLATION:

			if(degree == 2)
				localInterpQuadratic(points, size);
			else
				localInterpCubic(points, size);
			break;
	};
}

snlCurve::snlCurve(snlPoint& startPoint, snlPoint& endPoint, snlPoint& centrePoint, int numSections)
{
	// Create curve as circular arc.
	// -----------------------------
	// startPoint:      Starting point of arc.
	// endPoint:        Ending point of arc.
	// centrePoint:     Centre point of arc.
	// numSections:    (optional) specify number of sections arc has.
	//
	// Notes:           Can not do an arc bigger than 180 degrees.

	snlVector arcStart = startPoint - centrePoint;
	snlVector arcEnd = endPoint - centrePoint;

	double arcAngle = arcStart.angle(arcEnd);

	if(numSections == -1)
		numSections = (int)(( arcAngle /(M_PI / 2.0)) + 1);

	double sectionAngle = arcAngle / (double) numSections;

	double stepAngle = sectionAngle / 2.0;

	double midPtWeight = cos(stepAngle);

	int numCtrlPts =(numSections - 1) * 2 + 3;

	snlCtrlPoint* ctrlPts = new snlCtrlPoint[numCtrlPts];

	ctrlPts[0] = startPoint;
	ctrlPts[numCtrlPts - 1] = endPoint;

	snlTransform rotation;

	snlVector normal;

	// Rotate points into place.

	normal.crossProduct(arcStart, arcEnd);

	rotation.rotate(stepAngle, centrePoint, normal);

	for(int ptNum = 1; ptNum < numCtrlPts - 1; ptNum ++)
	{
		ctrlPts[ptNum] = ctrlPts[ptNum - 1];

		rotation.multiply(ctrlPts[ptNum]);
	}

	// Set mid point lengths and weights.

	int index = 1;

	// Calculate length of vector to add to mid point. midPtWeight is the cosine of the step angle.
	double midPtVectLength =(arcStart.length() / midPtWeight) - arcStart.length();

	for(int sectNum = 0; sectNum < numSections; sectNum ++)
	{
		snlVector midPtVect = ctrlPts[index] - centrePoint;

		midPtVect.length(midPtVectLength);

		ctrlPts[index] += midPtVect;

		ctrlPts[index].multiplyWeight(midPtWeight);

		index += 2;
	}

	// Generate control point net.

	ctrlPtNet = new snlCtrlPointNetCurve(ctrlPts, numCtrlPts);

	// Generate knot vector. Degree 2.

	int numKnots = numCtrlPts + 3;

	knot* knots = new knot[numKnots];

	// End clamps.

	for(index = 0; index < 3; index ++)
	{
		knots[index] = 0.0;
		knots[numKnots - index - 1] = 1.0;
	}

	// Internal knots.

	index = 3;

	for(int sectNum = 1; sectNum < numSections; sectNum ++)
	{
		knot knotVal =(1.0 / (double) numSections) * (double) sectNum;

		knots[index ++] = knotVal;
		knots[index ++] = knotVal;
	}

	// Generate knot vector.

	knotVect = new snlKnotVector(knots, numKnots, 2);

	deg = 2;
}

snlCurve& snlCurve::operator=(const snlCurve& curveToCopy)
{
	if(this != &curveToCopy)
	{
		if(ctrlPtNet) delete ctrlPtNet;
		if(knotVect) delete knotVect;

		ctrlPtNet = new snlCtrlPointNetCurve(*(curveToCopy.ctrlPtNet));

		knotVect = new snlKnotVector(*(curveToCopy.knotVect));

		deg = curveToCopy.deg;
	}

	return *this;
}

snlCtrlPointNetCurve& snlCurve::controlPointNet()
{
	// Return reference to control point network object for curve.
	// -----------------------------------------------------------

	return *ctrlPtNet;
}

snlPoint snlCurve::evalHmg(knot param) const
{
	// Evaluate Non Rational Homogeneous Curve Point.
	// ----------------------------------------------
	// param: Parameter to evaluate.
	//
	// Returns: Homogeneous point on curve.

	snlPoint rPoint;  // Return point.

	unsigned span = knotVect -> findSpan(param);

	// Evaluate basis functions.
	// Just fix the size of the evaluated basis function array to the max it can be.
	basis bVals[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];
	knotVect -> evalBasis(param, bVals);

	unsigned baseIndex = span - (unsigned) deg;

	// Get control point array.
	const snlCtrlPoint* ctrlPts = ctrlPtNet -> getCtrlPts();

	rPoint.zero();  // Set everything to zero.

	for(int index = 0; index <= deg; index ++)
	{
		// This is more efficient than the operator based equivalent becaues there is not construction of temp vectors.
		rPoint.multAdd(bVals[index], ctrlPts[baseIndex + index]);
	}

	return rPoint;
}

snlPoint snlCurve::eval(knot param) const
{
	// Evaluate rational non-homogeneous curve point.
	// ----------------------------------------------
	// param:       Parameter to evaluate.
	//
	// Returns:     Non-homogeneous point on curve.

	snlPoint retPoint = evalHmg(param);
	retPoint.project();

	return retPoint;
}

snlPoint* snlCurve::evalDerivsHmg(knot param, unsigned deriv) const
{
	// Evaluate Non Rational Homogeneous Surface Derivatives.
	// ------------------------------------------------------
	// param:       Parameter to evaluate at.
	// deriv        Derivative order to evaluate.
	//
	// Returns:     Array of snlPoint[deriv + 1]. Calling function
	//              must delete[] this array.

	snlPoint* retPnts = new snlPoint[deriv + 1];

	// Find spans
	unsigned span = knotVect -> findSpan(param);

	// Evaluate basis functions.
	// Just fix the size of the evaluated basis function array to the max it can be.
	basis bVals[SNL_KNOT_VECTOR_MAX_DEG * SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];
	knotVect -> evalBasisDeriv(param, deriv, bVals);

	unsigned baseIndex = span - (unsigned) deg;

	// Get control point array.
	const snlCtrlPoint* ctrlPts = ctrlPtNet -> getCtrlPts();

	for(unsigned derivIndex = 0; derivIndex <= deriv; derivIndex ++)
	{
		retPnts[derivIndex].zero();  // Set everything to zero.

		for(int index = 0; index <= deg; index ++)
		{
			retPnts[derivIndex].multAdd(bVals[index + derivIndex *(deg + 1)], ctrlPts[baseIndex + index]);
		}
	}

	return retPnts;
}

snlPoint* snlCurve::evalDerivs(knot param, unsigned deriv) const
{
	// Evaluate Rational Non-Homogeneous Surface Derivatives
	// -----------------------------------------------------
	// param:       Parameter to evaluate at.
	// deriv:       Derivative order to evaluate.
	//
	// Returns:     Array of snlPoint[deriv + 1]. Calling function must
	//              delete[] array.

	// Get homogeneous derivatives.
	snlPoint* derivPts = evalDerivsHmg(param, deriv);

	// Array for returning points in.
	snlPoint* evalPts = new snlPoint[deriv + 1];

	evalPts[0] = derivPts[0];  // First point in array is not a derivative.

	evalPts[0].project();

	double w0 = derivPts[0].w();

	snlVector sum;

	for(unsigned derivIndex = 1; derivIndex <= deriv; derivIndex ++)
	{
		sum.zero();

		for (unsigned index = 1; index <= derivIndex; index ++)
		{
			sum.multAdd((binCoefs::binCoefArray[derivIndex][index] * derivPts[index].w()), evalPts[derivIndex - index]);
		}

		evalPts[derivIndex] = derivPts[derivIndex] - sum;
		evalPts[derivIndex] *= w0;
		evalPts[derivIndex].w(0.0);  // Point is actually a 3D vector so w must always be 0.
	}

	delete[] derivPts;

	return evalPts;
}

snlVector snlCurve::velocity(knot param)
{
	// Velocity(first derivative) of curve.
	// ---------------------------------------
	// param:    Parameter to get velocity at.

	snlPoint* derivPoints = evalDerivs(param, 1);

	snlVector retVect(derivPoints[1]);

	delete[] derivPoints;

	return retVect;
}

void snlCurve::insertKnot(knot iParam, bool reallocate)
{
	// Insert a knot into knot vector and calculate new control points.
	// ---------------------------------------------------------------
	// iParam:      Parameter value to insert.
	// reallocate:  Reallocate memory for control points.
	//
	// Notes:       ctrlPts MUST have an additional point space allocated at the end of
	//              each line in the array for the chosen direction.

	unsigned        count, index;
	snlCtrlPoint    pt1, pt2;

	if(reallocate)
		ctrlPtNet -> grow();

	// Span new knot belongs to.
	unsigned span = knotVect -> findSpan(iParam);

	// Pre calculate alphas.
	double* alpha = new double[deg];

	for(count = 0; count < (unsigned) deg; count ++)
	{
		index = span - deg + 1 + count;
		alpha[count]  =(iParam -(knotVect -> getKnotVal(index)))
						   /(knotVect -> getKnotVal(index + deg) - knotVect -> getKnotVal(index));
	}

	// Build temp array to store new array of control points in.
	snlCtrlPoint* tmpPts = new snlCtrlPoint[deg];

	// Get pointer to control points.
	snlCtrlPoint* ctrlPts = ctrlPtNet -> getCtrlPtsPtr();

	// Calculate new control points.
	for(count = 0; count < (unsigned) deg; count ++)
	{
		index = span - deg + 1 + count;

		// Get first and second ctrl points to process with

		pt1 = ctrlPts[index];
		pt2 = ctrlPts[index - 1];

		pt1 *= alpha[count];
		tmpPts[count] = pt1;
		pt2 *=(1.0 - alpha[count]);
		tmpPts[count] += pt2;
	}

	// Place new points into array.

	// Copy non-altered control points forward one position at the end of the array.
	for(count =(ctrlPtNet -> size()) - 1; count > span; count --)
		ctrlPts[count] = ctrlPts[count - 1];

	// Copy new control points into array.
	for(count = 0; count < (unsigned) deg; count ++)
	{
		index = span - deg + 1 + count;
		ctrlPts[index] = tmpPts[count];
	}

	// Insert new knot into knot vector
	knotVect -> insertKnot(iParam);

	delete[] tmpPts;
	delete[] alpha;
}

void snlCurve::insertKnots(knot iParam, int numToInsert, bool reallocate)
{
	// Insert multiple knots.
	// ----------------------
	// iParam:        Parameter to insert.
	// numToInsert:   Number of knots to insert.
	// reallocate:    Reallocate memory for control points.

	for(int index = 0; index < numToInsert; index ++)
		insertKnot(iParam, reallocate);
}

double snlCurve::removeKnots(int numKnots, unsigned removalIndex, double tolerance)
{
	// Remove multiple knots from index.
	// ---------------------------------
	// numKnots:            Number of knots to remove.
	// removalIndex:        Index to remove knot from.
	// tolerance:           Maximum error allowed before knot removal aborted.
	//                      No tolerance if equals 0.
	//
	// Returns:             Tolerance achieved during knot removal whether successful or not.
	//
	// Notes:               Only removes multiples of the same parameter value initially at removal index.

	if(numKnots < 1) return 0.0;

	double maxTol = 0.0;

	double param = knotVect -> getKnotVal(removalIndex);

	int multi = knotVect -> findMultiplicity(removalIndex);

	int numToRemove = numKnots > multi ? multi : numKnots;

	for(int count = 0; count < numToRemove; count ++)
	{
		double tol = removeKnot(removalIndex, tolerance);

		if(tol > maxTol) maxTol = tol;

		removalIndex = knotVect -> findSpan(param);
	}

	return maxTol;
}

double snlCurve::removeKnot(unsigned removalIndex, double tolerance)
{
	// Remove knot from curve.
	// -----------------------
	// removalIndex:        Index to remove knot from.
	// tolerance:           Maximum error allowed before knot removal aborted.
	//                      No tolerance if equals 0.
	//
	// Returns:             Tolerance achieved during knot removal whether successful or not.

	knot rParam = knotVect -> getKnotVal(removalIndex);

	// Span knot to be removed belongs to. This will always adjust the removal index to
	// point to a non-zero span. ie Multiplicity check.
	unsigned rSpan = knotVect -> findSpan(rParam);

	// Find multiplicity of knot at index.
	unsigned multi = knotVect -> findMultiplicity(rSpan);

	// Calculate the number of equations.
	unsigned numEqns = deg - multi + 1;

	// Pre calculate alphas.
	double* alpha = knotVect -> calcRemovalAlphas(rSpan);

	// Build temp array to store new set of control points in.
	// First and last points are not new.
	snlCtrlPoint* tmpPts = new snlCtrlPoint[numEqns + 1];

	// Get control point array and calculate starting point for processing new points within it.

	const snlCtrlPoint* ctrlPts = ctrlPtNet -> getCtrlPts();

	unsigned ctrlPtIndex = rSpan - deg - 1;

	// Seed temp array.

	tmpPts[0] = ctrlPts[ctrlPtIndex ++];

	// Generate new points.

	for(unsigned index = 1; index <= numEqns; index ++)
	{
		tmpPts[index] = ctrlPts[ctrlPtIndex ++];
		tmpPts[index].multSub((1.0 - alpha[index - 1]), tmpPts[index - 1]);
		tmpPts[index] /= alpha[index - 1];
	}

	// If error is under tolerance then copy new points into control point array.

	double error = (tmpPts[numEqns] - ctrlPts[ctrlPtIndex]).length();

	if(error <= tolerance || tolerance == 0.0)
	{
		// Use original curve control point instead of newly created one in last of equations.
		tmpPts[numEqns] = ctrlPts[ctrlPtIndex];

		// Replace points in control point array.

		ctrlPtNet -> replacePoints(tmpPts + 1, numEqns, rSpan - deg, numEqns + 1);

		ctrlPtNet -> truncatePointSpace(1);

		knotVect -> removeKnot(rSpan);
	}

	// Clean up.

	delete[] alpha;
	delete[] tmpPts;

	return error;

}

void snlCurve::refine(double tolerance)
{
	// Refine control point net until tolerance is achieved.
	// -----------------------------------------------------

	if(deg <= 1) return;  // Degree 1 curves are straight lines.

	bool tolOk = false;

	while(! tolOk)
	{
		tolOk = true;

		for(int index = 0; (unsigned) index <(ctrlPtNet -> size()) - deg; index ++)
		{
			// Test for flatness

			double  flatness;

			flatness = ctrlPtNet -> calcFlatness(index, deg + 1);

			if(flatness > tolerance)
			{
				// Insert knot into surface. Half way between existing knots.

				int insertIndex = index + deg;

				knot insertParam =(( knotVect -> getKnotVal(insertIndex + 1)
									   - knotVect -> getKnotVal(insertIndex)) / 2)
									   + knotVect -> getKnotVal(insertIndex);

				insertKnot(insertParam, true);

				tolOk = false;

				index ++;  // If this is not done then nothing converges if the curvature is too great.
			}
		}
	}
}

double snlCurve::maxParam() const
{
	// Return maximum parameter value for curve.
	// -----------------------------------------

	return knotVect -> max();
}

double snlCurve::minParam() const
{
	// Return minimum parameter value for curve.
	// -----------------------------------------

	return knotVect -> min();
}

double snlCurve::param(unsigned index) const
{
	// Return parameter at specified knot index.
	// -----------------------------------------

	return knotVect -> getKnotVal(index);
}

int snlCurve::size()
{
	return ctrlPtNet -> size();
}

void snlCurve::truncate(knot param, bool keepLast, bool reparam)
{
	// Truncate curve.
	// ---------------
	// param:       Parameter to truncate at.
	// keepLast:    Keep last part of curve instead of first part.
	// reparam:     Reparameterise curve to original pararemeter boundaries.

	knot paramStart = knotVect -> min();
	knot paramEnd = knotVect -> max();

	if(param == paramStart || param == paramEnd) return;

	insertPartition(param);

	// Remove unwanted control points.

	unsigned span = knotVect -> findSpan(param);

	if(keepLast)
		span -= deg;
	else
		span = knotVect -> previousSpan(span);

	ctrlPtNet -> truncate(span, keepLast);

	// Remove unwanted knots.
	knotVect -> truncate(param, keepLast);

	// Reparameterise if required.

	if(reparam)
		reparameterise(paramStart, paramEnd);
}

void snlCurve::insertPartition(knot param)
{
	// Insert partition into curve.
	// ----------------------------
	// param:    Parameter to insert partition into.
	//
	// Notes:    Function basically makes sure degree knots are present
	//           at supplied parameter.

	int numToInsert = deg - knotVect -> findMultiplicity(param);

	for(int index = 0; index < numToInsert; index ++)
		insertKnot(param, true);
}

void snlCurve::reparameterise(knot startKnot, knot endKnot)
{
	// Do a linear Reparameterise on curve.
	// ------------------------------------
	// startKnot:    New starting knot value of knot vector.
	// endKnot:      Ending knot value of knot vector.
	//
	// Notes:        Linear reparameterisations don't effect control points.

	knotVect -> reparameterise(startKnot, endKnot);
}

void snlCurve::reverseEvalDirection()
{
	// Reverse curves parametric evaluation direction.
	// -----------------------------------------------

	// Reverse knot vector.

	knotVect -> reverse();

	// Reverse control points.

	ctrlPtNet -> reverse();
}

void snlCurve::globalInterpClosedLoop(int type, snlPoint* points, unsigned size, int degree, knot** retParams)
{
	// Global interpolation as closed loop.
	// ------------------------------------
	// type:      Type of global interpolation from SNL_FITTING_TYPES.
	// points:    Points to interpolate between.
	// size:      Number of points.
	// degree:    Resultant curve should be this degree.
	// retParams: Pointer to pointer that points to array of parameters that correspond to given points.
	//
	// Notes:     Array returned via retParams should be deleted by calling function.

	// Make sure first and last points aren't the same.

	if(points[0] == points[size - 1]) size --;

	// Create new points array with overlap to interpolate with.

	int newSize = size + degree * 2 + 1;

	snlPoint* newPoints = new snlPoint[newSize];

	// Starting overlap.

	int newIndex = 0;

	for(unsigned index = size - degree; index < size; index ++)
		newPoints[newIndex ++] = points[index];

	// Middle points.

	for(unsigned index = 0; index < size; index ++)
		newPoints[newIndex ++] = points[index];

	// Ending join and overlap.

	for(int index = 0; index < degree + 1; index ++)
		newPoints[newIndex ++] = points[index];

	// Pass new points to global interpolation function.

	knot* newRetParams;

	globalInterp(type, newPoints, newSize, degree, &newRetParams);

	// Truncate curve.

	knot paramStart = knotVect -> min();
	knot paramEnd = knotVect -> max();

	knot newParamStart = newRetParams[degree];
	knot newParamEnd = newRetParams[newSize - degree - 1];

	truncate(newParamStart, true, false);
	truncate(newParamEnd, false, false);

	// Reparameterise the curve and create new point return parameters.

	reparameterise( paramStart, paramEnd);

	if(retParams)
	{
		knot* params = new knot[size + 1];

		int paramIndex = degree;  // Discard overlap params.

		knot oldSpan = newParamEnd - newParamStart;
		knot newSpan = paramEnd - paramStart;

		// There is an additional parameter now because of the start and end points coinciding.

		for(unsigned index = 0; index <= size; index ++)
		{
			params[index] =(((newRetParams[paramIndex] - newParamStart) / oldSpan) * newSpan) + paramStart;
			paramIndex ++;
		}

		*retParams = params;
	}

	// Clean up.

	delete[] newRetParams;
	delete[] newPoints;

}

void snlCurve::globalInterp(int type, snlPoint* points, unsigned size, int degree, knot** retParams)
{
	// Global interpolation.
	// ---------------------
	// type:      Type of global interpolation from SNL_FITTING_TYPES.
	// points:    Points to interpolate between.
	// size:      Number of points.
	// degree:    Resultant curve should be this degree.
	// retParams: Pointer to pointer that points to array of parameters that correspond to given points.
	//
	// Notes:     Array returned via retParams should be deleted by calling function.

	if(knotVect) delete knotVect;
	if(ctrlPtNet) delete ctrlPtNet;

	deg = degree;

	// Generate parameters.

	knot* params = new knot[size];

	knot totalChordLength = 0.0;

	snlVector chord;

	// Intermediate step. Calculate (square root of) chord length.

	for(unsigned index = 1; index < size; index ++)
	{
		chord.diff(points[index - 1], points[index]);

		knot chordLength;

		if(type == SNL_GLOBAL_INTERP_CENTRIFUGAL)
			chordLength = sqrt(chord.length());
		else
			chordLength = chord.length();

		totalChordLength += chordLength;

		params[index] = chordLength;
	}

	// Calculate final parameter values.

	params[0] = 0.0;
	params[size - 1] = 1.0;

	for(unsigned index = 1; index < size - 1; index ++)
		params[index] = params[index - 1] + params[index] / totalChordLength;

	// Generate knot vector.

	knot* knots = new knot[size + degree + 1];

	unsigned index;

	// Start clamp.
	for(index = 0; index <= (unsigned) degree; index ++)
		knots[index] = 0.0;

	// End clamp.
	for(index = size; index < size + degree + 1; index ++)
		knots[index] = 1.0;

	// Internal knots.
	for(index = 1; index < size - degree; index ++)
	{
		knot sum = 0.0;

		for(unsigned paramIndex = index; paramIndex < index + degree; paramIndex ++)
			sum += params[paramIndex];

		knots[index + degree] = sum / degree;
	}

	knotVect = new snlKnotVector(knots, size + degree + 1, degree);

	// Setup and solve linear equations.

	// Generate coefficient array.

	unsigned arraySize = size * size;

	double* coeffs = new double[arraySize];

	// Zero everything to begin with.

	for(index = 0; index < arraySize; index ++)
		coeffs[index] = 0.0;

	// First and last rows just relfect clamps.

	coeffs[0] = 1.0;
	coeffs[arraySize - 1] = 1.0;

	// Fill middle rows with basis function values.

	basis basisVals[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];

	for(unsigned row = 1; row < size - 1; row ++)
	{
		knotVect -> evalBasis(params[row], basisVals);

		unsigned span = knotVect -> findSpan(params[row]);

		unsigned rowStartIndex = row * size;

		index = 0;

		for(unsigned col = span - degree; col <= span; col ++)
			coeffs[rowStartIndex + col] = basisVals[index ++];
	}

	// Generate right hand sides.

	double* rhSides = new double[size * 4];  // x, y, z and w.

	for(unsigned row = 0; row < size; row ++)
	{
		index = row * 4;

		rhSides[index] = points[row].components[0];  // x.
		rhSides[index + 1] = points[row].components[1];  // y.
		rhSides[index + 2] = points[row].components[2];  // z.
		rhSides[index + 3] = points[row].components[3];  // w.
	}

	// Pass data to linear solver.

	snlSquareLinear solver(size, 4, coeffs, rhSides);  // Constructor also calls snlSquareLinear::solve().

	// Copy solved points into new control points.

	snlCtrlPoint* ctrlPts = new snlCtrlPoint[size];

	for(unsigned row = 0; row < size; row ++)
	{
		index = row * 4;

		ctrlPts[row].components[0] = rhSides[index];
		ctrlPts[row].components[1] = rhSides[index + 1];
		ctrlPts[row].components[2] = rhSides[index + 2];
		ctrlPts[row].components[3] = rhSides[index + 3];
	}

	// Create control point net.

	ctrlPtNet = new snlCtrlPointNetCurve(ctrlPts, size, false);

	// Return parameters of given points if needed.

	if(retParams)
		*retParams = params;
	else
		delete[] params;
}

snlCtrlPoint* snlCurve::genGlobalInterpPoints(snlPoint* points, unsigned size, knot* params, snlKnotVector* knots)
{
	// Generate control points for global interpolation.
	// -------------------------------------------------
	// points:      Array of points to interpolate between.
	// size:        Size of points array.
	// params:      Parameters that correspond to points array.
	// knots:       Knot vector to use.
	//
	// Returns:     Array of control points the same size as the given "points" array.

	// *** Setup and solve linear equations ***

	// Generate coefficient array.

	unsigned arraySize = size * size;

	double* coeffs = new double[arraySize];

	// Zero everything to begin with.

	unsigned index;

	for(index = 0; index < arraySize; index ++)
		coeffs[index] = 0.0;

	// First and last rows just relfect clamps.

	coeffs[0] = 1.0;
	coeffs[arraySize - 1] = 1.0;

	// Fill middle rows with basis function values.

	basis basisVals[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];

	int deg = knots -> getDegree();

	for(unsigned row = 1; row < size - 1; row ++)
	{
		knots -> evalBasis(params[row], basisVals);

		unsigned span = knots -> findSpan(params[row]);

		unsigned rowStartIndex = row * size;

		index = 0;

		for(unsigned col = span - deg; col <= span; col ++)
			coeffs[rowStartIndex + col] = basisVals[index ++];
	}

	// Generate right hand sides.

	double* rhSides = new double[size * 4];  // x, y, z and w.

	for(unsigned row = 0; row < size; row ++)
	{
		index = row * 4;

		rhSides[index] = points[row].components[0];  // x.
		rhSides[index + 1] = points[row].components[1];  // y.
		rhSides[index + 2] = points[row].components[2];  // z.
		rhSides[index + 3] = points[row].components[3];  // w.
	}

	// Pass data to linear solver.

	snlSquareLinear solver(size, 4, coeffs, rhSides);  // Constructor also calls snlSquareLinear::solve().

	// Copy solved points into new control points.

	snlCtrlPoint* ctrlPts = new snlCtrlPoint[size];

	for(unsigned row = 0; row < size; row ++)
	{
		index = row * 4;

		ctrlPts[row].components[0] = rhSides[index];
		ctrlPts[row].components[1] = rhSides[index + 1];
		ctrlPts[row].components[2] = rhSides[index + 2];
		ctrlPts[row].components[3] = rhSides[index + 3];
	}

	return ctrlPts;
}

void snlCurve::localInterpQuadratic(snlPoint* points, unsigned size)
{
	// Local quadratic interpolation.
	// ------------------------------


}

void snlCurve::localInterpCubic(snlPoint* points, unsigned size)
{
	// Local cubic interpolation.
	// --------------------------


}

void snlCurve::synchronise(snlCurve& curve)
{
	//! Synchronise this curves knot vector to given curve.
	//  ---------------------------------------------------
	//! @param curve Curve to snync to.
	//!
	//! @par Notes: Knots are only ever added NOT removed.
	//!             So if curve has less multiplicity at a particular span index
	//!             then true synchronisation will not occur and the caller
	//!             should call the synchronise function on curve with this object
	//!             as it's argument.

	if(curve.deg != deg) return;  // Curve to sync to must have same degree.

	unsigned numSpans = curve.knotVect -> numSpans();

	if(numSpans < 2) return;

	// If the degree is the same then the first span will always have the same multiplicity for both curves.

	unsigned spanIndex = curve.knotVect -> firstSpan();
	spanIndex = curve.knotVect -> nextSpan(spanIndex);

	for(unsigned index = 1; index < numSpans; index ++)
	{
		knot param = curve.knotVect -> getKnotVal(spanIndex);

		int multi = curve.knotVect -> findMultiplicity(spanIndex);

		unsigned insertSpan = knotVect -> findSpan(param);  // Where param would be inserted in this object.

		// If knot already exists in this curve then reduce the number of knots inserted.

		if(knotVect -> getKnotVal(insertSpan) == param) multi -= knotVect -> findMultiplicity(insertSpan);

		if(multi > 0) insertKnots(param, multi, true);

		// Get next span.

		spanIndex = curve.knotVect -> nextSpan(spanIndex);
	}
}

void snlCurve::makeCompatible(snlCurve* curve)
{
	// Make this curve and given curve compatible.
	// -------------------------------------------
	// curve:   Curve to make this curve compatible with.

	// Make sure the degree of each curve is the same.

	if(deg > curve -> deg)
		curve -> elevateDegree(deg - curve -> deg);

	if(curve -> deg > deg)
		elevateDegree(curve -> deg - deg);

	// Make parametric ranges equal.

	knot thisMin = knotVect -> min();
	knot thisMax = knotVect -> max();

	knot compMin =(curve -> knotVect) -> min();  // Given curve.
	knot compMax =(curve -> knotVect) -> max();

	if(thisMin != compMin || thisMax != compMax)
	{
		// Reparameterise both curves.

		knot newMin = thisMin > compMin ? compMin : thisMin;
		knot newMax = thisMax > compMax ? thisMax : compMax;

		reparameterise(newMin, newMax);
		curve -> reparameterise(newMin, newMax);
	}

	// Synchronise the knot vectors of each curve.

	synchronise(*curve);
	curve -> synchronise(*this);
}

unsigned snlCurve::createBezierSegments(int** retNumKnotsAdded)
{
	// Create Bezier Segments over entire curve.
	// -----------------------------------------
	// retNumKnotsAdded:    Pointer to pointer to return array with number of inserted knots in it.
	//                      Caller must delete this array. First index in array corresponds to second knot span.
	//
	// Returns:             Number of elements in returned array.

	// Find number of knots to be inserted and reallocate curve memory in one pass.

	// Find number of non-zero spans.
	unsigned numSpans = knotVect -> numSpans();

	// Find first spans knot index.
	unsigned knotIndex = knotVect -> firstSpan();

	// Resize control points array just once for all knot insertions.

	unsigned span = knotIndex;
	unsigned extraKnots = 0;

	int* addedKnots = new int[numSpans - 1];  // First index corresponds to second span as no knots are added to first span.

	// Find amount to resize by.
	for(unsigned spanIndex = 1; spanIndex < numSpans; spanIndex ++)
	{
		span = knotVect -> nextSpan(span);

		addedKnots[spanIndex - 1] = deg -(knotVect -> findMultiplicity(span));

		extraKnots += addedKnots[spanIndex - 1];
	}

	// Append extra control point space to end of current control points.
	ctrlPtNet -> appendPointSpace(extraKnots);

	// Find knot index of second knot span.
	span = knotVect -> nextSpan(knotIndex);

	for(unsigned spanIndex = 0; spanIndex < numSpans - 1; spanIndex ++)
	{
		// Increase multiplicity of span to degree.

		knot insertParam = knotVect -> getKnotVal(span);

		insertKnots(insertParam, addedKnots[spanIndex], false);

		// Re-adjust current span index to account for inserted knots.
		span = knotVect -> nextSpan(span + addedKnots[spanIndex]);
	}

	*retNumKnotsAdded = addedKnots;

	return numSpans - 1;
}

void snlCurve::elevateBezierSegmentPointsDegree(int origDegree, int byDegree, const snlCtrlPoint* origPoints,
												  snlCtrlPoint* newPoints)
{
	// Calculate new control points for Bezier segment that is being degree elevated.
	// ------------------------------------------------------------------------------
	// origDegree:  Original degree of segment.
	// byDegree:    Number of degrees segment is being elevated by.
	// origPoints:  Original points of non-elevated segment. Expected size is origDegree + 1.
	// newPoints:   Array to store new points in. Must be size origDegree + byDegree + 1.

	int numNewPts = origDegree + byDegree;

	for(int index = 0; index <= numNewPts; index ++)
		newPoints[index].zero();

	for(int index = 0; index <= numNewPts; index ++)
	{
		int sumStart =(index - byDegree) > 0 ?(index - byDegree) : 0;
		int sumEnd = origDegree < index ? origDegree : index;

		for(int index2 = sumStart; index2 <= sumEnd; index2 ++)
		{
			double multiplier =((double) binCoefs::binCoefArray[origDegree][index2]) *
				((double) binCoefs::binCoefArray[byDegree][index - index2]) /
				((double) binCoefs::binCoefArray[origDegree + byDegree][index]);

			newPoints[index].multAdd(multiplier, origPoints[index2]);
		}
	}
}

void snlCurve::elevateDegree(int byDegree)
{
	// Elevate degree of curve
	// -----------------------
	// byDegree:    Number of degrees to elevate by.

	// Convert curve into Bezier segments.

	int* addedKnots;

	unsigned numSegments = createBezierSegments(&addedKnots);

	numSegments ++;  // Number returned is array size which is one less than number of segments.

	// Grow control point net.

	ctrlPtNet -> appendPointSpace(numSegments * byDegree);

	// Elevate degree of Bezier segments.

	int newSegmentSize = deg + byDegree + 1;

	snlCtrlPoint* tmpPts = new snlCtrlPoint[newSegmentSize];

	const snlCtrlPoint* ctrlPts = ctrlPtNet -> getCtrlPts();

	int ptsIndex = 0;

	unsigned spanIndex = deg * 2;

	// Generate new points per segment.

	for(unsigned segment = 0; segment < numSegments; segment ++)
	{
		elevateBezierSegmentPointsDegree(deg, byDegree, ctrlPts + ptsIndex, tmpPts);

		// Replace points in control point array. First and last points are not altered.
		ctrlPtNet -> replacePoints(tmpPts + 1, newSegmentSize - 2, ptsIndex + 1, deg - 1);

		ptsIndex += deg + byDegree;

		// Add knots to knot vector.
		knotVect -> increaseMultiplicity(spanIndex, byDegree);

		spanIndex += deg + byDegree;
	}

	// Make sure start clamp is of degree + 1 multiplicity.

	knotVect -> increaseMultiplicity(deg, byDegree);

	// Increase degree indicator variables

	deg += byDegree;

	knotVect -> setDegree(deg);

	// Remove number of knots that were added during knot insertion.

	spanIndex = knotVect -> firstSpan();

	spanIndex += deg;

	for(unsigned segment = 0; segment < numSegments - 1; segment ++)
	{
		removeKnots(addedKnots[segment], spanIndex, 0.0);

		spanIndex += deg - addedKnots[segment];
	}

	// Clean up.

	delete[] addedKnots;
	delete[] tmpPts;
}

void snlCurve::appendCurve(snlCurve* curveToAppend, bool copy)
{
	// Append a curve to this curve.
	// -----------------------------
	// curveToAppend:    Curve to append to this curve.
	// copy:             Copy curveToAppend and don't modify it in any way.
	//
	// Notes:            This curve may have it's degree elevated to accomodate new curve.


	// Copy curve if required.

	snlCurve* curve;

	if(copy)
		curve = new snlCurve(*curveToAppend);
	else
		curve = curveToAppend;

	// Elevate degree if needed.

	if(deg > curve -> degree())
		curve -> elevateDegree(deg - curve -> degree());
	else if(deg < curve -> degree())
		elevateDegree(curve -> degree() - deg);

	// Re-parameterise curve to append so that it's starting knot val is the same as this curves knot end val.

	double min = knotVect -> min();
	double max = knotVect -> max();

	curve -> reparameterise(max, max + 1.0);

	// Join control points and knot vectors together.

	ctrlPtNet -> appendPoints(( curve -> ctrlPtNet) -> getCtrlPts() + 1,(curve -> ctrlPtNet) -> getNumPts() - 1);

	knotVect -> join(curve -> knotVect);

	reparameterise(min, max);

	// Clean up.

	if(copy) delete curve;
}

void snlCurve::print()
{
	// Print curve to standard out.
	// ----------------------------

	ctrlPtNet -> print();
	knotVect -> print();
}

void snlCurve::vertexNet(snlVertexNet* vNet, double tolerance, bool parametric)
{
	// Return approximation to curve.
	// --------------------------------
	// tolerance:    Tolerance to approximate to.
	// parametric:   Do a parametric analysis as opposed to knot refinement.
	// vNet:         Vertex net to fill with data.

	if(parametric)
	{
		vertexNetParam(vNet, tolerance);
		return;
	}

	const snlCtrlPoint*   ctrlPts;
	int                   size;

	snlCurve* tmpCurve = 0;

	if(tolerance > 0.0)
	{
		tmpCurve = new snlCurve(*this);
		tmpCurve -> refine(tolerance);
		ctrlPts =(tmpCurve -> ctrlPtNet) -> getCtrlPts();
		size =(tmpCurve -> ctrlPtNet) -> size();
	}
	else
	{
		ctrlPts = ctrlPtNet -> getCtrlPts();
		size = ctrlPtNet -> size();
	}

	vNet -> vertexNet(ctrlPts, size);

	if(tmpCurve) delete tmpCurve;
}

void snlCurve::vertexNetParam(snlVertexNet* vNet, double tolerance)
{
	// Generate vertex net based on evaluation of curve.
	// -------------------------------------------------
	// vNet:        Vertex network to load with points.
	// tolerance:   Tolerance to actual surface.

	int              size;

	snlPoint* pts = 0;

	if(tolerance <= 0.0)
	{
		size = ctrlPtNet -> size();

		pts = new snlPoint[size];

		double minP = minParam();

		double paramStep =(maxParam() - minP) / (double)(size - 1);

		for(int index = 0; index < size; index ++)
		{
			double param = minP + paramStep * (double) index;
			pts[index] = eval(param);
		}
	}
	else
	{
		// !@#$ Not complete!!
	}

	if(pts) vNet -> vertexNet(pts, size);

	if(pts) delete[] pts;
}

const snlKnotVector& snlCurve::knotVector()
{
	return *knotVect;
}

int snlCurve::degree()
{
	return deg;
}

