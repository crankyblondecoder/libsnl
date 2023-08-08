// *** General NURBS Curve ***

#include "snlCurve.h"
#include "snlUtil.h"
#include "snlSquareLinear.h"

snlCurve::~snlCurve()
{
	if(_ctrlPtNet) delete _ctrlPtNet;
	if(_knotVect) delete _knotVect;
}

snlCurve::snlCurve()
{
	_deg = 0;

	_ctrlPtNet = 0;
	_knotVect = 0;
}

snlCurve::snlCurve(const snlCurve& copyFrom)
{
	_ctrlPtNet = new snlCtrlPointNetCurve(*(copyFrom._ctrlPtNet));

	_knotVect = new snlKnotVector(*(copyFrom._knotVect));

	_deg = copyFrom._deg;
}

snlCurve& snlCurve::operator=(const snlCurve& curveToCopy)
{
	if(this != &curveToCopy)
	{
		if(_ctrlPtNet) delete _ctrlPtNet;
		if(_knotVect) delete _knotVect;

		_ctrlPtNet = new snlCtrlPointNetCurve(*(curveToCopy._ctrlPtNet));

		_knotVect = new snlKnotVector(*(curveToCopy._knotVect));

		_deg = curveToCopy._deg;
	}

	return *this;
}

snlCurve::snlCurve(int degree, unsigned size, snlPoint& start, snlPoint& end)
{
	_deg = degree;

	_ctrlPtNet = new snlCtrlPointNetCurve(size, start, end);

	_knotVect = new snlKnotVector(0.0, 1.0, size + degree + 1, degree);
}

snlCurve::snlCurve(int degree, unsigned size, snlCtrlPoint* points, knot* knots)
{
	_deg = degree;

	_ctrlPtNet = new snlCtrlPointNetCurve(points, size, false);

	if(knots)
		_knotVect = new snlKnotVector(knots, size + degree + 1, degree);
	else
		_knotVect = new snlKnotVector(0.0, 1.0, size + degree + 1, degree);
}

snlCurve::snlCurve(unsigned size, snlCtrlPoint* points, snlKnotVector* knotVector)
{
	_deg = knotVector -> getDegree();

	_knotVect = knotVector;

	_ctrlPtNet = new snlCtrlPointNetCurve(points, size, false);
}

snlCurve::snlCurve(snlPoint* points, unsigned size, int fittingType, int degree, bool closedLoop, knot** retParams)
{
	_ctrlPtNet = 0;
	_knotVect = 0;

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

	_ctrlPtNet = new snlCtrlPointNetCurve(ctrlPts, numCtrlPts);

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

	_knotVect = new snlKnotVector(knots, numKnots, 2);

	_deg = 2;
}

snlCtrlPointNetCurve& snlCurve::controlPointNet()
{
	return *_ctrlPtNet;
}

snlPoint snlCurve::evalHmg(knot param) const
{
	snlPoint rPoint;  // Return point.

	unsigned span = _knotVect -> findSpan(param);

	// Evaluate basis functions.
	// Just fix the size of the evaluated basis function array to the max it can be.
	basis bVals[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];
	_knotVect -> evalBasis(param, bVals);

	unsigned baseIndex = span - (unsigned) _deg;

	// Get control point array.
	const snlCtrlPoint* ctrlPts = _ctrlPtNet -> getCtrlPts();

	rPoint.zero();  // Set everything to zero.

	for(int index = 0; index <= _deg; index ++)
	{
		// This is more efficient than the operator based equivalent becaues there is not construction of temp vectors.
		rPoint.multAdd(bVals[index], ctrlPts[baseIndex + index]);
	}

	return rPoint;
}

snlPoint snlCurve::eval(knot param) const
{
	snlPoint retPoint = evalHmg(param);
	retPoint.project();

	return retPoint;
}

snlPoint* snlCurve::evalDerivsHmg(knot param, unsigned deriv) const
{
	snlPoint* retPnts = new snlPoint[deriv + 1];

	// Find spans
	unsigned span = _knotVect -> findSpan(param);

	// Evaluate basis functions.
	// Just fix the size of the evaluated basis function array to the max it can be.
	basis bVals[SNL_KNOT_VECTOR_MAX_DEG * SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];
	_knotVect -> evalBasisDeriv(param, deriv, bVals);

	unsigned baseIndex = span - (unsigned) _deg;

	// Get control point array.
	const snlCtrlPoint* ctrlPts = _ctrlPtNet -> getCtrlPts();

	for(unsigned derivIndex = 0; derivIndex <= deriv; derivIndex ++)
	{
		retPnts[derivIndex].zero();  // Set everything to zero.

		for(int index = 0; index <= _deg; index ++)
		{
			retPnts[derivIndex].multAdd(bVals[index + derivIndex *(_deg + 1)], ctrlPts[baseIndex + index]);
		}
	}

	return retPnts;
}

snlPoint* snlCurve::evalDerivs(knot param, unsigned deriv) const
{
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
	snlPoint* derivPoints = evalDerivs(param, 1);

	snlVector retVect(derivPoints[1]);

	delete[] derivPoints;

	return retVect;
}

void snlCurve::insertKnot(knot iParam, bool reallocate)
{
	unsigned        count, index;
	snlCtrlPoint    pt1, pt2;

	if(reallocate)
		_ctrlPtNet -> grow();

	// Span new knot belongs to.
	unsigned span = _knotVect -> findSpan(iParam);

	// Pre calculate alphas.
	double* alpha = new double[_deg];

	for(count = 0; count < (unsigned) _deg; count ++)
	{
		index = span - _deg + 1 + count;
		alpha[count]  =(iParam -(_knotVect -> getKnotVal(index)))
						   /(_knotVect -> getKnotVal(index + _deg) - _knotVect -> getKnotVal(index));
	}

	// Build temp array to store new array of control points in.
	snlCtrlPoint* tmpPts = new snlCtrlPoint[_deg];

	// Get pointer to control points.
	snlCtrlPoint* ctrlPts = _ctrlPtNet -> getCtrlPtsPtr();

	// Calculate new control points.
	for(count = 0; count < (unsigned) _deg; count ++)
	{
		index = span - _deg + 1 + count;

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
	for(count =(_ctrlPtNet -> size()) - 1; count > span; count --)
		ctrlPts[count] = ctrlPts[count - 1];

	// Copy new control points into array.
	for(count = 0; count < (unsigned) _deg; count ++)
	{
		index = span - _deg + 1 + count;
		ctrlPts[index] = tmpPts[count];
	}

	// Insert new knot into knot vector
	_knotVect -> insertKnot(iParam);

	delete[] tmpPts;
	delete[] alpha;
}

void snlCurve::insertKnots(knot iParam, int numToInsert, bool reallocate)
{
	for(int index = 0; index < numToInsert; index ++)
		insertKnot(iParam, reallocate);
}

double snlCurve::removeKnots(int numKnots, unsigned removalIndex, double tolerance)
{
	if(numKnots < 1) return 0.0;

	double maxTol = 0.0;

	double param = _knotVect -> getKnotVal(removalIndex);

	int multi = _knotVect -> findMultiplicity(removalIndex);

	int numToRemove = numKnots > multi ? multi : numKnots;

	for(int count = 0; count < numToRemove; count ++)
	{
		double tol = removeKnot(removalIndex, tolerance);

		if(tol > maxTol) maxTol = tol;

		removalIndex = _knotVect -> findSpan(param);
	}

	return maxTol;
}

double snlCurve::removeKnot(unsigned removalIndex, double tolerance)
{
	knot rParam = _knotVect -> getKnotVal(removalIndex);

	// Span knot to be removed belongs to. This will always adjust the removal index to
	// point to a non-zero span. ie Multiplicity check.
	unsigned rSpan = _knotVect -> findSpan(rParam);

	// Find multiplicity of knot at index.
	unsigned multi = _knotVect -> findMultiplicity(rSpan);

	// Calculate the number of equations.
	unsigned numEqns = _deg - multi + 1;

	// Pre calculate alphas.
	double* alpha = _knotVect -> calcRemovalAlphas(rSpan);

	// Build temp array to store new set of control points in.
	// First and last points are not new.
	snlCtrlPoint* tmpPts = new snlCtrlPoint[numEqns + 1];

	// Get control point array and calculate starting point for processing new points within it.

	const snlCtrlPoint* ctrlPts = _ctrlPtNet -> getCtrlPts();

	unsigned ctrlPtIndex = rSpan - _deg - 1;

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

		_ctrlPtNet -> replacePoints(tmpPts + 1, numEqns, rSpan - _deg, numEqns + 1);

		_ctrlPtNet -> truncatePointSpace(1);

		_knotVect -> removeKnot(rSpan);
	}

	// Clean up.

	delete[] alpha;
	delete[] tmpPts;

	return error;

}

void snlCurve::refine(double tolerance)
{
	if(_deg <= 1) return;  // Degree 1 curves are straight lines.

	bool tolOk = false;

	while(! tolOk)
	{
		tolOk = true;

		for(int index = 0; (unsigned) index <(_ctrlPtNet -> size()) - _deg; index ++)
		{
			// Test for flatness

			double  flatness;

			flatness = _ctrlPtNet -> calcFlatness(index, _deg + 1);

			if(flatness > tolerance)
			{
				// Insert knot into surface. Half way between existing knots.

				int insertIndex = index + _deg;

				knot insertParam =(( _knotVect -> getKnotVal(insertIndex + 1)
									   - _knotVect -> getKnotVal(insertIndex)) / 2)
									   + _knotVect -> getKnotVal(insertIndex);

				insertKnot(insertParam, true);

				tolOk = false;

				index ++;  // If this is not done then nothing converges if the curvature is too great.
			}
		}
	}
}

double snlCurve::maxParam() const
{
	return _knotVect -> max();
}

double snlCurve::minParam() const
{
	return _knotVect -> min();
}

double snlCurve::param(unsigned index) const
{
	return _knotVect -> getKnotVal(index);
}

int snlCurve::size()
{
	return _ctrlPtNet -> size();
}

void snlCurve::truncate(knot param, bool keepLast, bool reparam)
{
	knot paramStart = _knotVect -> min();
	knot paramEnd = _knotVect -> max();

	if(param == paramStart || param == paramEnd) return;

	insertPartition(param);

	// Remove unwanted control points.

	unsigned span = _knotVect -> findSpan(param);

	if(keepLast)
		span -= _deg;
	else
		span = _knotVect -> previousSpan(span);

	_ctrlPtNet -> truncate(span, keepLast);

	// Remove unwanted knots.
	_knotVect -> truncate(param, keepLast);

	// Reparameterise if required.

	if(reparam)
		reparameterise(paramStart, paramEnd);
}

void snlCurve::insertPartition(knot param)
{
	int numToInsert = _deg - _knotVect -> findMultiplicity(param);

	for(int index = 0; index < numToInsert; index ++)
		insertKnot(param, true);
}

void snlCurve::reparameterise(knot startKnot, knot endKnot)
{
	_knotVect -> reparameterise(startKnot, endKnot);
}

void snlCurve::reverseEvalDirection()
{
	// Reverse knot vector.

	_knotVect -> reverse();

	// Reverse control points.

	_ctrlPtNet -> reverse();
}

void snlCurve::globalInterpClosedLoop(int type, snlPoint* points, unsigned size, int degree, knot** retParams)
{
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

	knot paramStart = _knotVect -> min();
	knot paramEnd = _knotVect -> max();

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
	if(_knotVect) delete _knotVect;
	if(_ctrlPtNet) delete _ctrlPtNet;

	_deg = degree;

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

	_knotVect = new snlKnotVector(knots, size + degree + 1, degree);

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

	basis basisVals[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];

	for(unsigned row = 1; row < size - 1; row ++)
	{
		_knotVect -> evalBasis(params[row], basisVals);

		unsigned span = _knotVect -> findSpan(params[row]);

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

	_ctrlPtNet = new snlCtrlPointNetCurve(ctrlPts, size, false);

	// Return parameters of given points if needed.

	if(retParams)
		*retParams = params;
	else
		delete[] params;
}

snlCtrlPoint* snlCurve::genGlobalInterpPoints(snlPoint* points, unsigned size, knot* params, snlKnotVector* knots)
{
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

	basis basisVals[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];

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
	// TODO ...
}

void snlCurve::localInterpCubic(snlPoint* points, unsigned size)
{
	// TODO ...
}

void snlCurve::synchronise(snlCurve& curve)
{
	if(curve._deg != _deg) return;  // Curve to sync to must have same degree.

	unsigned numSpans = curve._knotVect -> numSpans();

	if(numSpans < 2) return;

	// If the degree is the same then the first span will always have the same multiplicity for both curves.

	unsigned spanIndex = curve._knotVect -> firstSpan();
	spanIndex = curve._knotVect -> nextSpan(spanIndex);

	for(unsigned index = 1; index < numSpans; index ++)
	{
		knot param = curve._knotVect -> getKnotVal(spanIndex);

		int multi = curve._knotVect -> findMultiplicity(spanIndex);

		unsigned insertSpan = _knotVect -> findSpan(param);  // Where param would be inserted in this object.

		// If knot already exists in this curve then reduce the number of knots inserted.

		if(_knotVect -> getKnotVal(insertSpan) == param) multi -= _knotVect -> findMultiplicity(insertSpan);

		if(multi > 0) insertKnots(param, multi, true);

		// Get next span.

		spanIndex = curve._knotVect -> nextSpan(spanIndex);
	}
}

void snlCurve::makeCompatible(snlCurve* curve)
{
	// Make sure the degree of each curve is the same.

	if(_deg > curve -> _deg)
		curve -> elevateDegree(_deg - curve -> _deg);

	if(curve -> _deg > _deg)
		elevateDegree(curve -> _deg - _deg);

	// Make parametric ranges equal.

	knot thisMin = _knotVect -> min();
	knot thisMax = _knotVect -> max();

	knot compMin =(curve -> _knotVect) -> min();  // Given curve.
	knot compMax =(curve -> _knotVect) -> max();

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
	// Find number of knots to be inserted and reallocate curve memory in one pass.

	// Find number of non-zero spans.
	unsigned numSpans = _knotVect -> numSpans();

	// Find first spans knot index.
	unsigned knotIndex = _knotVect -> firstSpan();

	// Resize control points array just once for all knot insertions.

	unsigned span = knotIndex;
	unsigned extraKnots = 0;

	int* addedKnots = new int[numSpans - 1];  // First index corresponds to second span as no knots are added to first span.

	// Find amount to resize by.
	for(unsigned spanIndex = 1; spanIndex < numSpans; spanIndex ++)
	{
		span = _knotVect -> nextSpan(span);

		addedKnots[spanIndex - 1] = _deg -(_knotVect -> findMultiplicity(span));

		extraKnots += addedKnots[spanIndex - 1];
	}

	// Append extra control point space to end of current control points.
	_ctrlPtNet -> appendPointSpace(extraKnots);

	// Find knot index of second knot span.
	span = _knotVect -> nextSpan(knotIndex);

	for(unsigned spanIndex = 0; spanIndex < numSpans - 1; spanIndex ++)
	{
		// Increase multiplicity of span to degree.

		knot insertParam = _knotVect -> getKnotVal(span);

		insertKnots(insertParam, addedKnots[spanIndex], false);

		// Re-adjust current span index to account for inserted knots.
		span = _knotVect -> nextSpan(span + addedKnots[spanIndex]);
	}

	*retNumKnotsAdded = addedKnots;

	return numSpans - 1;
}

void snlCurve::elevateBezierSegmentPointsDegree(int origDegree, int byDegree, const snlCtrlPoint* origPoints,
	snlCtrlPoint* newPoints)
{
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
	// Convert curve into Bezier segments.

	int* addedKnots;

	unsigned numSegments = createBezierSegments(&addedKnots);

	numSegments ++;  // Number returned is array size which is one less than number of segments.

	// Grow control point net.

	_ctrlPtNet -> appendPointSpace(numSegments * byDegree);

	// Elevate degree of Bezier segments.

	int newSegmentSize = _deg + byDegree + 1;

	snlCtrlPoint* tmpPts = new snlCtrlPoint[newSegmentSize];

	const snlCtrlPoint* ctrlPts = _ctrlPtNet -> getCtrlPts();

	int ptsIndex = 0;

	unsigned spanIndex = _deg * 2;

	// Generate new points per segment.

	for(unsigned segment = 0; segment < numSegments; segment ++)
	{
		elevateBezierSegmentPointsDegree(_deg, byDegree, ctrlPts + ptsIndex, tmpPts);

		// Replace points in control point array. First and last points are not altered.
		_ctrlPtNet -> replacePoints(tmpPts + 1, newSegmentSize - 2, ptsIndex + 1, _deg - 1);

		ptsIndex += _deg + byDegree;

		// Add knots to knot vector.
		_knotVect -> increaseMultiplicity(spanIndex, byDegree);

		spanIndex += _deg + byDegree;
	}

	// Make sure start clamp is of degree + 1 multiplicity.

	_knotVect -> increaseMultiplicity(_deg, byDegree);

	// Increase degree indicator variables

	_deg += byDegree;

	_knotVect -> setDegree(_deg);

	// Remove number of knots that were added during knot insertion.

	spanIndex = _knotVect -> firstSpan();

	spanIndex += _deg;

	for(unsigned segment = 0; segment < numSegments - 1; segment ++)
	{
		removeKnots(addedKnots[segment], spanIndex, 0.0);

		spanIndex += _deg - addedKnots[segment];
	}

	// Clean up.

	delete[] addedKnots;
	delete[] tmpPts;
}

void snlCurve::appendCurve(snlCurve* curveToAppend, bool copy)
{
	// Copy curve if required.

	snlCurve* curve;

	if(copy)
		curve = new snlCurve(*curveToAppend);
	else
		curve = curveToAppend;

	// Elevate degree if needed.

	if(_deg > curve -> degree())
		curve -> elevateDegree(_deg - curve -> degree());
	else if(_deg < curve -> degree())
		elevateDegree(curve -> degree() - _deg);

	// Re-parameterise curve to append so that it's starting knot val is the same as this curves knot end val.

	double min = _knotVect -> min();
	double max = _knotVect -> max();

	curve -> reparameterise(max, max + 1.0);

	// Join control points and knot vectors together.

	_ctrlPtNet -> appendPoints(( curve -> _ctrlPtNet) -> getCtrlPts() + 1,(curve -> _ctrlPtNet) -> getNumPts() - 1);

	_knotVect -> join(curve -> _knotVect);

	reparameterise(min, max);

	// Clean up.

	if(copy) delete curve;
}

void snlCurve::print()
{
	_ctrlPtNet -> print();
	_knotVect -> print();
}

void snlCurve::vertexNet(snlVertexNet* vNet, double tolerance, bool parametric)
{
	if(parametric)
	{
		_vertexNetParam(vNet, tolerance);
		return;
	}

	const snlCtrlPoint*   ctrlPts;
	int                   size;

	snlCurve* tmpCurve = 0;

	if(tolerance > 0.0)
	{
		tmpCurve = new snlCurve(*this);
		tmpCurve -> refine(tolerance);
		ctrlPts =(tmpCurve -> _ctrlPtNet) -> getCtrlPts();
		size =(tmpCurve -> _ctrlPtNet) -> size();
	}
	else
	{
		ctrlPts = _ctrlPtNet -> getCtrlPts();
		size = _ctrlPtNet -> size();
	}

	vNet -> vertexNet(ctrlPts, size);

	if(tmpCurve) delete tmpCurve;
}

void snlCurve::_vertexNetParam(snlVertexNet* vNet, double tolerance)
{
	int size;

	snlPoint* pts = 0;

	if(tolerance <= 0.0)
	{
		size = _ctrlPtNet -> size();

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
	return *_knotVect;
}

int snlCurve::degree()
{
	return _deg;
}

