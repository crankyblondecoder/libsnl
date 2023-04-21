#include "snlCtrlPointNetCurve.h"

snlCtrlPointNetCurve::snlCtrlPointNetCurve(snlCtrlPoint* cPtArray, unsigned size, bool copy)
{
	// Control Points for a curve - Constructor
	// ----------------------------------------
	// cPtArray:    Array of points to copy.
	// size_u:      Size of array to copy.
	// copy:        Make a copy of cPtArray.

	ctrlPtArraySize = size;

	if(copy)
	{
		// Copy points into object.
		ctrlPts = new snlCtrlPoint [ ctrlPtArraySize ];

		for(unsigned count = 0; count < ctrlPtArraySize; count ++)
			ctrlPts [ count ] = cPtArray [ count ];
	}
	else
		ctrlPts = cPtArray;
}

snlCtrlPointNetCurve::snlCtrlPointNetCurve(unsigned size, snlPoint& start, snlPoint& end)
{
	// Construct a control point network.
	// ----------------------------------
	// size:    Number of control points.
	// start:   Starting point of curve.
	// end:     Ending point of curve.

	if(! size)
	{
		ctrlPts = 0;
		return;
	}

	ctrlPtArraySize = size;

	ctrlPts = new snlCtrlPoint [ ctrlPtArraySize ];

	snlVector lineVect = end - start;

	lineVect *=(1.0 /(size - 1));

	for(unsigned index = 0; index < size; index ++)
	{
		snlVector offsetVect = lineVect * ((double) index);
		snlPoint newPoint = start + offsetVect;
		ctrlPts [ index ] = newPoint;  // Conversion to snlCtrlPoint
	}
}

snlCtrlPointNetCurve::~snlCtrlPointNetCurve()
{
}

unsigned snlCtrlPointNetCurve::size() const
{
	return ctrlPtArraySize;
}

snlCtrlPoint* snlCtrlPointNetCurve::grow()
{
	// Increase the control point net's size.
	// --------------------------------------

	if(! ctrlPts) return 0;

	snlCtrlPoint* newPts = new snlCtrlPoint [ ctrlPtArraySize + 1 ];

	// Copy points into new array.
	for(unsigned index = 0; index < ctrlPtArraySize; index ++)
		newPts [ index ] = ctrlPts [ index ];

	// Delete old array and point to new one.
	delete[] ctrlPts;

	ctrlPts = newPts;

	ctrlPtArraySize ++;

	return ctrlPts;
}

snlCtrlPoint* snlCtrlPointNetCurve::shrink()
{
	// Decrease the control point net's size.
	// --------------------------------------

	if(! ctrlPts) return 0;

	snlCtrlPoint* newPts = new snlCtrlPoint [ ctrlPtArraySize - 1 ];

	// Copy points into new array.
	for(unsigned index = 0; index <(ctrlPtArraySize - 1); index ++)
		newPts [ index ] = ctrlPts [ index ];

	// Delete old array and point to new one.

	delete[] ctrlPts;

	ctrlPts = newPts;

	ctrlPtArraySize --;

	return ctrlPts;
}

double snlCtrlPointNetCurve::calcFlatness(int index, int numPoints) const
{
	// Calculate flatness of a series of points.
	// -----------------------------------------
	// index:        Index to get points from.
	// numPoints:    Number of points to evaluate.

	if((unsigned)(index + numPoints) > ctrlPtArraySize) return 0;

	snlPoint** testPoints  = new snlPoint* [ numPoints ];

	for(int count = 0; count < numPoints; count ++)
	{
		testPoints [ count ] = ctrlPts + index + count;
	}

	double flatness = snlCtrlPointNet::calcFlatness(testPoints, numPoints);

	delete[] testPoints;

	return flatness;
}

void snlCtrlPointNetCurve::truncate(int atIndex, bool keepLast)
{
	// Truncate control point array.
	// -----------------------------
	// index:    Array index to truncate from.
	// keepLast: Keep last part of array and truncate first part.

	snlCtrlPoint* newCtrlPts;

	unsigned newSize;

	if(keepLast)
	{
		newSize = ctrlPtArraySize - atIndex;

		newCtrlPts = new snlCtrlPoint [ newSize ];

		for(unsigned index = 0; index < newSize; index ++)
			newCtrlPts [ index ] = ctrlPts [ atIndex + index ];
	}
	else
	{
		newSize = atIndex + 1;

		newCtrlPts = new snlCtrlPoint [ newSize ];

		for(unsigned index = 0; index < newSize; index ++)
			newCtrlPts [ index ] = ctrlPts [ index ];
	}

	ctrlPtArraySize = newSize;

	delete[] ctrlPts;

	ctrlPts = newCtrlPts;
}

void snlCtrlPointNetCurve::reverse()
{
	// Reverse order of control point net array.
	// -----------------------------------------

	unsigned midPoint = ctrlPtArraySize / 2;

	unsigned swapIndex = ctrlPtArraySize - 1;

	snlCtrlPoint trans;

	for(unsigned index = 0; index < midPoint; index ++)
	{
		trans = ctrlPts [ index ];  // Transfer value.
		ctrlPts [ index ] = ctrlPts [ swapIndex ];
		ctrlPts [ swapIndex -- ] = trans;
	}
}

int snlCtrlPointNetCurve::getCnctPts(unsigned index, snlCtrlPoint* retPts)
{
	// Return connected control points.
	// --------------------------------
	// index:    Index of point to find connections of.
	// retPts:   Array to return points in.
	//
	// Returns:  Number of connected points found.

	if(index >= ctrlPtArraySize) return 0;

	int retIndex = 0;

	if(index > 0)
	{
		retPts [ retIndex ] = ctrlPts [ index - 1 ];
		retIndex ++;
	}

	if(index <(ctrlPtArraySize - 1))
	{
		retPts [ retIndex ] = ctrlPts [ index + 1 ];
		retIndex ++;
	}

	return retIndex;
}

int snlCtrlPointNetCurve::maxConnections() const
{
	return 2;
}

