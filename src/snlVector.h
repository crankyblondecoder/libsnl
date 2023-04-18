#ifndef SNL_VECTOR_H
#define SNL_VECTOR_H

#include <iostream>
#include <cmath>

using namespace std;

/**
 * Vector in R⁴.
 */
class snlVector
{
	public:

		/**
		 * New vector which is zeroed.
		 */
		snlVector();

		/** New vector as deep copy from other vector. */
		snlVector(const snlVector& copyFrom);

		/**
		 * New vector given three or more components.
		 */
		snlVector(double x, double y, double z, double w = 0.0);

		/**
		 * Set the components of this as the vector which if added to the first vector gives the second vector.
		 * @param vec1 The first vector.
		 * @param vec2 The second vector.
		 */
		void diff(const snlVector& v1, const snlVector& v2);

		/**
		 * Calculate cross product of v1 X v2 in R³ and store in this vector.
		 * @note The w coordinate is ignored from the given vectors and set to zero in the result.
		 */
		void crossProduct(snlVector& v1, snlVector& v2);

		/**
		 * Calculate dot product between this vector and the given vector.
		 */
		double dot(snlVector& vect);

		/**
		 * Return length of this vector squared.
		 */
		double lengthSqrd();  // Length of vector squared.

		/**
		 * Calculate the length of this vector.
		 */
		double length();

		/**
		 * Set the length of this vector.
		 */
		void length(double length);

		/**
		 * Calculate absolute cosine of the angle between this and the given vector.
		 * ie. The absolute value of the cosine of the angle between the vectors.
		 */
		double calcAbsCos(snlVector& vect);

		/**
		 * Calculate angle between this vector and the given vector.
		 * @returns Value between 0 and PI in radians.
		 */
		double angle(snlVector& vect);

		/**
		 * Normalise this vector. ie Converts it to a unit vector.
		 */
		void normalise();

		/**
		 * Caclulate projection distance from tip of given vector to this vector.
		 * @param fromVector Vector to project tip of.
		 * @note Vectors are assumed to originate from the same point. This is standard projection as defined in linear
		 *       algebra text books.
		 */
		double projectDist(snlVector& fromVector);

		/**
		 * Project this vector onto another vector.
		 * @param ontoVector Vector to project onto.
		 * @returns New vector that is the result of the projection.
		 */
		snlVector project(snlVector& ontoVector);

// HERE
		void projectXZ();  // Project onto the X-Z plane.
		void projectXY();  // Project onto the X-Y plane.
		void projectYZ();  // Project onto the Y-Z plane.

		snlVector operator * (double);  // Return vector multiplied by a scalar.
		snlVector operator + (snlVector& vect);
		snlVector operator - (snlVector& vect);

		void operator += (snlVector& vect);
		void operator -= (snlVector& vect);
		void operator *= (double);  // Multiply this vector by a scalar.

		bool operator == (snlVector& compare);

		double x();
		double y();
		double z();
		double w();

		void x(double val);
		void y(double val);
		void z(double val);
		void w(double val);

		void components(double x, double y, double z, double w);
		void components(double x, double y, double z);

		void components(double* x, double* y, double* z, double* w);
		void components(double* x, double* y, double* z);

		void zero();  // Zero the vector.

		bool isZero();

		void print();

		double elements[4];
};

// Inline Functions
// ----------------

inline void snlVector::diff(const snlVector& v1, const snlVector& v2)
{
	elements[0] = v2.elements[0] - v1.elements[0];
	elements[1] = v2.elements[1] - v1.elements[1];
	elements[2] = v2.elements[2] - v1.elements[2];
	elements[3] = v2.elements[3] - v1.elements[3];
}

inline void snlVector::crossProduct (snlVector& v1, snlVector& v2)
{
	elements[0] = (v1.elements[1] * v2.elements[2]) - (v1.elements[2] * v2.elements[1]);
	elements[1] = (v1.elements[2] * v2.elements[0]) - (v1.elements[0] * v2.elements[2]);
	elements[2] = (v1.elements[0] * v2.elements[1]) - (v1.elements[1] * v2.elements[0]);
	elements[3] = 0.0;
}

inline double snlVector::dot (snlVector& vect)
{
	return (elements[0] * vect.elements[0] +
			elements[1] * vect.elements[1] +
			elements[2] * vect.elements[2] +
			elements[3] * vect.elements[3]);
}

inline double snlVector::lengthSqrd()
{
	return elements[0] * elements[0] +
		   elements[1] * elements[1] +
		   elements[2] * elements[2] +
		   elements[3] * elements[3];
}

#endif
