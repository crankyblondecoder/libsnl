#ifndef SNLPOINT_H
#define SNLPOINT_H

class snlVector;

#include "snlVector.h"

/**
 * Point in homogeneous coordinates (R⁴) that can be projected to R³.
 */
class snlPoint
{
	public:

		/** New point initialised to origin. */
		snlPoint();
		/** New point as deep copy from other point. */
		snlPoint(const snlPoint& copyFrom);
		/** New point with given coordinates. */
		snlPoint(double x, double y, double z, double w = 1.0);

		/** Set all components of this point.*/
		void setComponents(double x, double y, double z, double w);
		/** Set just the cartesian components of this point. W is unaffected. */
		void setComponents(double x, double y, double z);

		/** Get the homogeneous components of this point. Populates the locations of the given pointers.*/
		void getComponents(double* x, double* y, double* z, double* w) const;
		/** Get just the cartesian components of this point. Populates the locations of the given pointers. */
		void getComponents(double* x, double* y, double* z) const;

		/** Get the points x coordinate. */
		double x() const;
		/** Get the points y coordinate. */
		double y() const;
		/** Get the points z coordinate. */
		double z() const;
		/** Get the points w coordinate. */
		double w() const;

		/** Set the points x coordinate. */
		void x (double);
		/** Set the points y coordinate. */
		void y (double);
		/** Set the points z coordinate. */
		void z (double);
		/** Set the points w coordinate. */
		void w (double);

		/** Multiply this points weight (w coordinate) by a multiplier. */
		void multiplyWeight(double multiplier);

		/**
		 * Project homogeneous point into R³.
		 * ie Divide all coordinate, including w, by w.
		 */
		void project();

		/** Set all coordinates to zero. */
		void null();

		/**
		 * Get whether all coordinates are zero.
		 */
		bool isNull();

		/**
		 * Set point to [0,0,0,1]
		 * Essentially moves point to origin in R^3.
		 */
		void zeroInR3();

		/** Return length of this point, treated as a 4D vector, squared. */
		double lengthSqrd() const;

		/**
		 * Return squared distance from this point to given point.
		 * @param toPoint Point to calculate distance to.
		 * @note This is a 3D distance only. ie Distance calcs are only done with normalised points.
		 */
		double distSqrd ( const snlPoint& toPoint ) const;

		/** Add a vector to this point and return a new point. */
		snlPoint operator + ( const snlVector& vect ) const;
		/** Subtract a vector from this point and return a new point. */
		snlPoint operator - ( const snlVector& vect ) const;
		/** Subtract a point from this point and return a new vector. */
		snlVector operator - ( const snlPoint& point ) const;
		/** Multiply this point by a scalar and return a new point. */
		snlPoint operator * ( double scalar ) const;
		/** Divide this point by a scalar and return a new point. */
		snlPoint operator / ( double scalar ) const;

		void operator = ( const snlPoint& copyFrom );
		void operator += ( const snlPoint& point );
		void operator += ( const snlVector& vect );
		void operator -= ( const snlPoint& point );
		void operator *= ( double scalar );
		void operator /= ( double scalar );

		bool operator == ( snlPoint& compare );

		void print() const;

		/** [x, y, z, w] */
		double elements[4];
};

#endif
