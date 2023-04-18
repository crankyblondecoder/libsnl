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

void snlVector::print()
{
	// Print vector contents to cout.
	// ------------------------------

	cout << "X: " << elements[0] << " Y: " << elements[1]
		 << " Z: " << elements[2] << " W: " << elements[3] << "\n";
}
