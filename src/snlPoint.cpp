#include "snlPoint.h"

snlPoint::snlPoint()
{
	elements[0] = 0;
	elements[1] = 0;
	elements[2] = 0;
	elements[3] = 1;
}

snlPoint::snlPoint(const snlPoint& copyFrom)
{
	for(unsigned index = 0; index < 4; index ++)
	{
		elements[index] = copyFrom.elements[index];
	}
}

snlPoint::snlPoint(double x, double y, double z, double w)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = w;
}

void snlPoint::setComponents(double x, double y, double z, double w)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = w;
}

void snlPoint::setComponents(double x, double y, double z)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
}

void snlPoint::getComponents(double* x, double* y, double* z, double* w) const
{
	*x = elements[0];
	*y = elements[1];
	*z = elements[2];
	*w = elements[3];
}

void snlPoint::getComponents(double* x, double* y, double* z) const
{
	*x = elements[0];
	*y = elements[1];
	*z = elements[2];
}

double snlPoint::x() const
{
	return elements[0];
}

double snlPoint::y() const
{
	return elements[1];
}

double snlPoint::z() const
{
	return elements[2];
}

double snlPoint::w() const
{
	return elements[3];
}

void snlPoint::multiplyWeight(double multiplier)
{
	elements[0] *= multiplier;
	elements[1] *= multiplier;
	elements[2] *= multiplier;
	elements[3] *= multiplier;
}

void snlPoint::project()
{
	if(elements[3] == 0.0) return;  // Stop divide by zero error.

	double w = elements[3];

	elements[0] /= w;
	elements[1] /= w;
	elements[2] /= w;
	elements[3] = 1.0;
}

void snlPoint::null()
{
	elements[0] = 0;
	elements[1] = 0;
	elements[2] = 0;
	elements[3] = 0;
}

bool snlPoint::isNull()
{
	// Just rely on w only ever being 0 if this point is null.
	return elements[3] == 0;
}

void snlPoint::zeroInR3()
{
	elements[0] = 0;
	elements[1] = 0;
	elements[2] = 0;
	elements[3] = 1;
}

snlPoint snlPoint::operator + (const snlVector& vect) const
{
	// This should be copied on return.
	snlPoint retPt;

	if(vect.homogeneous)
	{
		retPt.elements[0] = elements[0] + vect.elements[0];
		retPt.elements[1] = elements[1] + vect.elements[1];
		retPt.elements[2] = elements[2] + vect.elements[2];
		retPt.elements[3] = elements[3] + vect.elements[3];
	}
	else
	{
		// Need to adjust the vector elements for the w coordinate. The vector is treated as though w = 1.
		retPt.elements[0] = elements[0] + vect.elements[0] * elements[3];
		retPt.elements[1] = elements[1] + vect.elements[1] * elements[3];
		retPt.elements[2] = elements[2] + vect.elements[2] * elements[3];
		retPt.elements[3] = elements[3];
	}

	return retPt;
}

snlPoint snlPoint::operator - (const snlVector& vect) const
{
	// This should be copied on return.
	snlPoint retPt;

	if ( vect.homogeneous )
	{
		retPt.elements[0] = elements[0] - vect.elements[0];
		retPt.elements[1] = elements[1] - vect.elements[1];
		retPt.elements[2] = elements[2] - vect.elements[2];
		retPt.elements[3] = elements[3] - vect.elements[3];
	}
	else
	{
		// Need to adjust the vector elements for the w coordinate. The vector is treated as though w = 1.
		retPt.elements[0] = elements[0] - vect.elements[0] * elements[3];
		retPt.elements[1] = elements[1] - vect.elements[1] * elements[3];
		retPt.elements[2] = elements[2] - vect.elements[2] * elements[3];
		retPt.elements[3] = elements[3] - vect.elements[3] * elements[3];
		retPt.elements[3] = elements[3];
	}

	return retPt;
}

snlVector snlPoint::operator - (const snlPoint& point) const
{
	// This should be copied on return.
	snlVector retVec;

	retVec.elements[0] = elements[0] - point.elements [0];
	retVec.elements[1] = elements[1] - point.elements [1];
	retVec.elements[2] = elements[2] - point.elements [2];
	retVec.elements[3] = elements[3] - point.elements [3];

	return retVec;
}

snlPoint snlPoint::operator * ( double scalar ) const
{
	// This should be copied on return.
	snlPoint retPt;

	retPt.elements[0] = elements[0] * scalar;
	retPt.elements[1] = elements[1] * scalar;
	retPt.elements[2] = elements[2] * scalar;
	retPt.elements[3] = elements[3] * scalar;

	return retPt;
}

snlPoint snlPoint::operator / (double scalar) const
{
	// This should be copied on return.
	snlPoint retPt;

	for ( int index = 0; index < 4; index ++ )
		retPt.elements [ index ] = elements [ index ] / scalar;

	return retPt;
}

void snlPoint::operator = ( const snlPoint& copyFrom )
{
	// Copy data from another point.
	// -----------------------------

	elements [ 0 ] = copyFrom.elements [ 0 ];
	elements [ 1 ] = copyFrom.elements [ 1 ];
	elements [ 2 ] = copyFrom.elements [ 2 ];
	elements[3] = copyFrom.elements[3];
}

void snlPoint::operator += ( const snlPoint& point )
{
	// Add a point to this point.
	// --------------------------

	elements [ 0 ] += point.elements [ 0 ];
	elements [ 1 ] += point.elements [ 1 ];
	elements [ 2 ] += point.elements [ 2 ];
	elements[3] += point.elements[3];
}

void snlPoint::operator += ( const snlVector& vect )
{
	// Add a vector to this point.
	// ---------------------------

	if ( vect.homogeneous )
	{
		elements [ 0 ] += vect.elements [ 0 ];
		elements [ 1 ] += vect.elements [ 1 ];
		elements [ 2 ] += vect.elements [ 2 ];
		elements[3] += vect.elements[3];
	}
	else
	{
		elements [ 0 ] += vect.elements [ 0 ] * elements[3];
		elements [ 1 ] += vect.elements [ 1 ] * elements[3];
		elements [ 2 ] += vect.elements [ 2 ] * elements[3];
	}
}

void snlPoint::operator *= ( double scalar )
{
	// Mulitply a scalar to this point.
	// --------------------------------

	for ( int index = 0; index < 4; index ++ )
		elements [ index ] *= scalar;
}

void snlPoint::operator /= ( double scalar )
{
	// Divide this point by a scalar.
	// ------------------------------

	for ( int index = 0; index < 4; index ++ )
		elements [ index ] /= scalar;
}

void snlPoint::x ( double val )
{
	// Set x component.
	// ----------------

	elements [ 0 ] = val;
}

void snlPoint::y ( double val )
{
	// Set y component.
	// ----------------

	elements [ 1 ] = val;
}

void snlPoint::z ( double val )
{
	// Set z component.
	// ----------------

	elements [ 2 ] = val;
}

void snlPoint::w ( double val )
{
	// Set w component.
	// ----------------

	elements[3] = val;
}

double snlPoint::lengthSqrd() const
{
	double sum = 0;

	for ( int index = 0; index < 4; index ++ )
		sum += elements [ index ] * elements [ index ];

	return sum;
}

double snlPoint::distSqrd(const snlPoint& toPoint) const
{
	double diff;
	double retVal = 0;

	if(elements[3] == 1.0 && toPoint.elements[3] == 1.0)
	{
		diff = elements[0] - toPoint.elements[0];
		retVal += diff * diff;
		diff = elements[1] - toPoint.elements[1];
		retVal += diff * diff;
		diff = elements[2] - toPoint.elements[2];
		retVal += diff * diff;
	}
	else
	{
		diff = (elements[0] / elements[3]) - (toPoint.elements[0] / toPoint.elements[3]);
		retVal += diff * diff;
		diff = (elements[1] / elements[3]) - (toPoint.elements[1] / toPoint.elements[3]);
		retVal += diff * diff;
		diff = (elements[2] / elements[3]) - (toPoint.elements[2] / toPoint.elements[3]);
		retVal += diff * diff;
	}

	return retVal;
}

bool snlPoint::operator == ( snlPoint& compare )
{
	// Return true if compare is eqivalent to this point.
	// --------------------------------------------------

	bool retVal = true;

	for ( int index = 0; index < 4; index ++ )
		if ( elements [ index ] != compare.elements [ index ] )
			retVal = false;

	return retVal;
}

void snlPoint::print() const
{
	// Print to std out.
	// -----------------

	cout << "( " << elements [ 0 ] << ", " <<  elements [ 1 ] << ", " << elements [ 2 ] << ", " << elements[3] << " )";
}

