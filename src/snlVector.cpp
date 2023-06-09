#include "snlVector.h"

snlVector::snlVector()
{
	components[0] = 0.0;
	components[1] = 0.0;
	components[2] = 0.0;
	components[3] = 0.0;
}

snlVector::snlVector(double x, double y, double z, double w)
{
	components[0] = x;
	components[1] = y;
	components[2] = z;
	components[3] = w;
}

double snlVector::length()
{
	return sqrt(lengthSqrd());
}

void snlVector::length(double len)
{
	double multiplier = len / length();

	components[0] *= multiplier;
	components[1] *= multiplier;
	components[2] *= multiplier;
	components[3] *= multiplier;
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

	components[0] /= len;
	components[1] /= len;
	components[2] /= len;
	components[3] /= len;
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

	cout << "X: " << components[0] << " Y: " << components[1]
		 << " Z: " << components[2] << " W: " << components[3] << "\n";
}
