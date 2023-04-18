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

		/**
		 * Project onto the X-Z plane.
		 */
		void projectXZ();

		/**
		 * Project onto the X-Y plane.
		 */
		void projectXY();

		/**
		 * Project onto the Y-Z plane.
		 */
		void projectYZ();

		/** Return new vector that is this vector multiplied by a scalar. */
		snlVector operator * (double);
		/** Return new vector that is this vector plus another vector. */
		snlVector operator + (snlVector& vect);
		/** Return new vector that is this vector minus another vector. */
		snlVector operator - (snlVector& vect);

		/** Add a vector to this vector. */
		void operator += (snlVector& vect);
		/** Subract a vector from this vector. */
		void operator -= (snlVector& vect);
		/** Multiply this vector by a scalar. */
		void operator *= (double);

		bool operator == (snlVector& compare);

		/** Get the x coordinate. */
		double x();
		/** Get the y coordinate. */
		double y();
		/** Get the z coordinate. */
		double z();
		/** Get the w coordinate. */
		double w();

		/** Set the x coordinate. */
		void x(double val);
		/** Set the y coordinate. */
		void y(double val);
		/** Set the z coordinate. */
		void z(double val);
		/** Set the w coordinate. */
		void w(double val);

		/** Set the vectors x, y, z and w components. */
		void components(double x, double y, double z, double w);
		/** Set the vectors x, y and z components. */
		void components(double x, double y, double z);

		/** Get the vectors x, y, z and w components by copying to given pointers. */
		void components(double* x, double* y, double* z, double* w);
		/** Get the vectors x, y and z components by copying to given pointers. */
		void components(double* x, double* y, double* z);

		/** Set this vector to the zero vector, ie [0,0,0,0] */
		void zero();

		/** Get whether this vector is the zero vector. */
		bool isZero();

		void print();

	protected:

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

inline void snlVector::projectXZ()
{
	elements[1] = 0.0;  // y = 0.
}

inline void snlVector::projectXY()
{
	elements[2] = 0.0;  // z = 0.
}

inline void snlVector::projectYZ()
{
	elements[0] = 0.0;
}

inline double snlVector::x()
{
	return elements[0];
}

inline double snlVector::y()
{
	return elements[1];
}

inline double snlVector::z()
{
	return elements[2];
}

inline double snlVector::w()
{
	return elements[3];
}

inline void snlVector::x(double val)
{
	elements[0] = val;
}

inline void snlVector::y(double val)
{
	elements[1] = val;
}

inline void snlVector::z(double val)
{
	elements[2] = val;
}

inline void snlVector::w(double val)
{
	elements[3] = val;
}

inline void snlVector::components(double x, double y, double z, double w)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = w;
}

inline void snlVector::components(double x, double y, double z)
{
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	elements[3] = 0.0;
}

inline void snlVector::components(double* x, double* y, double* z, double* w)
{
	*x = elements[0];
	*y = elements[1];
	*z = elements[2];
	*w = elements[3];
}

inline void snlVector::components(double* x, double* y, double* z)
{
	*x = elements[0];
	*y = elements[1];
	*z = elements[2];
}

inline void snlVector::zero()
{
	elements[0] = 0.0;
	elements[1] = 0.0;
	elements[2] = 0.0;
	elements[3] = 0.0;
}

inline bool snlVector::isZero()
{

	return !(
		elements[0] ||
		elements[1] ||
		elements[2] ||
		elements[3]
	);
}

inline snlVector snlVector::operator * (double scalar)
{
	return snlVector(
		elements[0] * scalar,
		elements[1] * scalar,
		elements[2] * scalar,
		elements[3] * scalar);
}

inline snlVector snlVector::operator + (snlVector& vect)
{
	return	snlVector(
		elements[0] + vect.elements[0],
		elements[1] + vect.elements[1],
		elements[2] + vect.elements[2],
		elements[3] + vect.elements[3]);
}

inline snlVector snlVector::operator - (snlVector& vect)
{
	return	snlVector(
		elements[0] - vect.elements[0],
		elements[1] - vect.elements[1],
		elements[2] - vect.elements[2],
		elements[3] - vect.elements[3]);
}

inline void snlVector::operator += (snlVector& vect)
{
	elements[0] += vect.elements[0];
	elements[1] += vect.elements[1];
	elements[2] += vect.elements[2];
	elements[3] += vect.elements[3];
}

inline void snlVector::operator -= (snlVector& vect)
{
	elements[0] -= vect.elements[0];
	elements[1] -= vect.elements[1];
	elements[2] -= vect.elements[2];
	elements[3] -= vect.elements[3];
}

inline void snlVector::operator *= (double scalar)
{
	elements[0] *= scalar;
	elements[1] *= scalar;
	elements[2] *= scalar;
	elements[3] *= scalar;
}

inline bool snlVector::operator == (snlVector& compare)
{
	return
		elements[0] == compare.elements[0] &&
		elements[1] == compare.elements[1] &&
		elements[2] == compare.elements[2] &&
		elements[3] == compare.elements[3];
}

#endif
