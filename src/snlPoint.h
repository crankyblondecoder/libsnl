#ifndef SNLPOINT_H
#define SNLPOINT_H

class snlVector;

#include "snlVector.h"

/**
 * Point in homogeneous coordinates (R⁴) that can be projected to R³.
 * Represented as a position vector.
 */
class snlPoint : public snlVector
{
	public:

		/** New point initialised to origin. */
		snlPoint();

		/** New point with given coordinates. */
		snlPoint(double x, double y, double z, double w = 1.0);

		/** Multiply this points weight (w coordinate) by a multiplier. */
		void multiplyWeight(double multiplier);

		/** Calculate the distance squared between this and another point. */
		double distSqrd(const snlPoint& point);

		/** Calculate the distance between this and another point. */
		double distance(const snlPoint& point);

		/**
		 * Project homogeneous point into R³.
		 * ie Divide all coordinate, including w, by w.
		 */
		void project();

		/**
		 * Set point to [0,0,0,1]
		 * Essentially moves point to origin in R³.
		 */
		void zeroInR3();

		/** Return new Point that is this Point plus a Vector. */
		snlPoint operator + (const snlVector& vect) const;

		/** Return new Point that is this Point minus a Vector. */
		snlPoint operator - (const snlVector& vect) const;
};

inline void snlPoint::multiplyWeight(double multiplier)
{
	elements[0] *= multiplier;
	elements[1] *= multiplier;
	elements[2] *= multiplier;
	elements[3] *= multiplier;
}

inline double snlPoint::distSqrd(const snlPoint& point)
{
	// This is a duplicate of combined vector operations but is necessary for speed.

	double retVal = 0.0;
	double scratch;

	scratch = point.elements[0] - elements[0];
	retVal += scratch * scratch;

	scratch = point.elements[1] - elements[1];
	retVal += scratch * scratch;

	scratch = point.elements[2] - elements[2];
	retVal += scratch * scratch;

	scratch = point.elements[3] - elements[3];
	retVal += scratch * scratch;

	return retVal;
}

inline double snlPoint::distance(const snlPoint& point)
{
	double dSqrd = distSqrd(point);
	return sqrt(dSqrd);
}

inline void snlPoint::project()
{
	if(elements[3] == 0.0) return;  // Stop divide by zero error.

	double w = elements[3];

	elements[0] /= w;
	elements[1] /= w;
	elements[2] /= w;
	elements[3] = 1.0;
}

inline void snlPoint::zeroInR3()
{
	elements[0] = 0;
	elements[1] = 0;
	elements[2] = 0;
	elements[3] = 1;
}

#endif
