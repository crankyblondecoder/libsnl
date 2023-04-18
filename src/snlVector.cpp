#include "snlVector.h"

snlVector::snlVector()
{
	elements[0] = 0.0;
	elements[1] = 0.0;
	elements[2] = 0.0;
	elements[3] = 0.0;
}

snlVector::snlVector(const snlVector& copyFrom)
{
	elements[0] = copyFrom.elements[0];
	elements[1] = copyFrom.elements[1];
	elements[2] = copyFrom.elements[2];
	elements[3] = copyFrom.elements[3];
}

snlVector::snlVector(double x, double y, double z, double w)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = w;
}

double snlVector::dot(snlPoint& pt)
{
	// Calculate dot product between this vector and pt.
	// -------------------------------------------------
	// Notes:   Treats pt as a vector.

	double dotProd = elements[0] * pt.elements[0] +
					 elements[1] * pt.elements[1] +
					 elements[2] * pt.elements[2];

	if(homogeneous)
		dotProd += elements[3] * pt.elements[3];

	return dotProd;
}


snlVector snlVector::operator *(double scalar)
{
	// Return vector multiplied by a scalar
	// ------------------------------------

	snlVector    retVect;

	for(int index = 0; index < 4; index ++)
		retVect.elements[index] = elements[index] * scalar;

	return retVect;
}

snlVector snlVector::operator +(snlVector& vect)
{
	// Add vect to this one.
	// ---------------------

	snlVector    retVect;

	for(int index = 0; index < 4; index ++)
		retVect.elements[index] = elements[index] + vect.elements[index];

	return retVect;
}

snlVector snlVector::operator -(snlVector& vect)
{
	// Subtract vect from this one.
	// ----------------------------

	snlVector    retVect;

	for(int index = 0; index < 4; index ++)
		retVect.elements[index] = elements[index] - vect.elements[index];

	return retVect;
}

void snlVector::operator +=(snlVector& vect)
{
	elements[0] += vect.elements[0];
	elements[1] += vect.elements[1];
	elements[2] += vect.elements[2];

	if(homogeneous)
		elements[3] += vect.elements[3];
}

void snlVector::operator -=(snlVector& vect)
{
	elements[0] -= vect.elements[0];
	elements[1] -= vect.elements[1];
	elements[2] -= vect.elements[2];

	if(homogeneous)
		elements[3] -= vect.elements[3];
}

void snlVector::operator *=(double scalar)
{
	// Multiply this vector by scalar.
	// -------------------------------

	for(int index = 0; index < 4; index ++)
		elements[index] *= scalar;
}

bool snlVector::operator ==(snlVector& compare)
{
	// Return true if compare is eqivalent to this vector.
	// ---------------------------------------------------

	bool retVal = true;

	if(homogeneous)
	{
		for(int index = 0; index < 4; index ++)
			if(elements[index] != compare.elements[index])
				retVal = false;
	}
	else
	{
		for(int index = 0; index < 3; index ++)
			if(elements[index] != compare.elements[index])
				retVal = false;
	}

	return retVal;
}

double snlVector::length()
{
	return sqrt(lengthSqrd());
}

void snlVector::length(double len)
{
	double multiplier = len / length();

	elements[0] *= multiplier;
	elements[1] *= multiplier;
	elements[2] *= multiplier;
	elements[3] *= multiplier;
}

double snlVector::calcAbsCos(snlVector& vect)
{
	return fabs(dot(vect) / (length() * vect.length()));
}

double snlVector::angle(snlVector& vect)
{
	double cos_angle = dot(vect) /(length() * vect.length());

	// Adjust for floating point error.

	if(cos_angle > 1.0)
	{
		cos_angle = 1.0;
	}
	else if(cos_angle < -1.0)
	{
		cos_angle = -1.0;
	}

	return acos(cos_angle);
}

void snlVector::normalise()
{
	double len = length();

	elements[0] /= len;
	elements[1] /= len;
	elements[2] /= len;
	elements[3] /= len;
}

double snlVector::projectDist(snlVector& fromVector)
{
	double dotP = dot(fromVector);

	return sqrt(fromVector.lengthSqrd() -(dotP * dotP / lengthSqrd()));
}

snlVector snlVector::project(snlVector& ontoVector)
{
	double newLength = dot(ontoVector) / ontoVector.length();

	snlVector retVector = ontoVector;

	retVector.length(newLength);

	return retVector;
}

void snlVector::projectXZ()
{
	// Project onto the X-Z plane.
	// ---------------------------

	elements[1] = 0.0;  // y = 0.
}

void snlVector::projectXY()
{
	// Project onto the X-Y plane.
	// ---------------------------

	elements[2] = 0.0;  // z = 0.
}

void snlVector::projectYZ()
{
	// Project onto the Y-Z plane.
	// ---------------------------

	elements[0] = 0.0;
}

double snlVector::x()
{
	return elements[0];
}

double snlVector::y()
{
	return elements[1];
}

double snlVector::z()
{
	return elements[2];
}

double snlVector::w()
{
	return elements[3];
}

void snlVector::x(double val)
{
	elements[0] = val;
}

void snlVector::y(double val)
{
	elements[1] = val;
}

void snlVector::z(double val)
{
	elements[2] = val;
}

void snlVector::w(double val)
{
	elements[3] = val;
}

void snlVector::calcNormal(snlPoint& pt1, snlPoint& pt2, snlPoint& pt3, snlPoint& pt4)
{
	// Calculate normal and set this vector to it.
	// -------------------------------------------
	// NOTE:        Only does 3D vector product.

	snlVector v1(pt1, pt2);
	snlVector v2(pt3, pt4);

	crossProduct(v1, v2);

	unitise();
}

void snlVector::components(double x, double y, double z, double w)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = w;
}

void snlVector::components(double x, double y, double z)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = 0.0;
}

void snlVector::components(double* x, double* y, double* z, double* w)
{
	*x = elements[0];
	*y = elements[1];
	*z = elements[2];
	*w = elements[3];
}

void snlVector::components(double* x, double* y, double* z)
{
	*x = elements[0];
	*y = elements[1];
	*z = elements[2];
}

void snlVector::zero()
{
	// Zero the vector.
	// ----------------

	elements[0] = 0.0;
	elements[1] = 0.0;
	elements[2] = 0.0;
	elements[3] = 0.0;
}

bool snlVector::isZero()
{
	if(homogeneous)
	{
		if(! elements[0] &&
			 ! elements[1] &&
			 ! elements[2] &&
			 ! elements[3])
			return true;
	}
	else
	{
		if(! elements[0] &&
			 ! elements[1] &&
			 ! elements[2])
		return true;
	}

	return false;
}

void snlVector::print()
{
	// Print vector contents to cout.
	// ------------------------------

	cout << "X: " << elements[0] << " Y: " << elements[1]
		 << " Z: " << elements[2] << " W: " << elements[3] << "\n";
}
