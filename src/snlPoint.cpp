#include "snlPoint.h"

snlPoint::snlPoint()
	: snlVector(0.0, 0.0, 0.0, 1.0)
{
}

snlPoint::snlPoint(double x, double y, double z, double w)
	: snlVector(x, y, z, w)
{
}

snlPoint snlPoint::operator + (snlVector& vect)
{
	snlPoint newPoint(*this);
	newPoint += vect;
	return newPoint;
}

snlPoint snlPoint::operator - (snlVector& vect)
{
	snlPoint newPoint(*this);
	newPoint -= vect;
	return newPoint;
}
