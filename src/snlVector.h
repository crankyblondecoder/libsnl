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

		/**
		 * New vector given three or more components.
		 */
		snlVector(double x, double y, double z, double w);

		/**
		 * Copy the elements from another vector into this.
		 */
		void assign(const snlVector& copyFrom);

		/**
		 * Set the components of this as the vector which if added to the first vector gives the second vector.
		 * ie vec2 - vec1.
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
		 * @note The calculated value takes into account the w component so unless these are zero the output of the function
		 *       might be unexpected.
		 */
		double calcAbsCos(snlVector& vect);

		/**
		 * Calculate angle between this vector and the given vector.
		 * @returns Value between 0 and PI in radians.
		 * @note The calculated value takes into account the w component so unless these are zero the output of the function
		 *       might be unexpected.
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
		snlVector operator * (double scalar) const;
		/** Return new vector that is the cross product of the given vector and this vector. */
		snlVector operator * (const snlVector& vect) const;
		/** Return new vector that is this vector plus another vector. */
		snlVector operator + (const snlVector& vect) const;
		/** Return new vector that is this vector minus another vector. */
		snlVector operator - (const snlVector& vect) const;

		/** Add a vector to this vector. */
		snlVector& operator += (const snlVector& vect);
		/** Subract a vector from this vector. */
		snlVector& operator -= (const snlVector& vect);
		/** Multiply this vector by a scalar. */
		snlVector& operator *= (double);
		/** Divide this vector by a scalar. */
		snlVector& operator /= (double);

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
		void setComponents(double x, double y, double z, double w);
		/** Set the vectors x, y and z components. */
		void setComponents(double x, double y, double z);

		/** Get the vectors x, y, z and w components by copying to given pointers. */
		void getComponents(double* x, double* y, double* z, double* w);
		/** Get the vectors x, y and z components by copying to given pointers. */
		void getComponents(double* x, double* y, double* z);

		/** Set this vector to the zero vector, ie [0,0,0,0] */
		void zero();

		/** Get whether this vector is the zero vector. */
		bool isZero();

		/**
		 * Multply the addend by a scalar, without modifying it, then add to this vector.
		 */
		void multAdd(double scalar, const snlVector& addend);

		/**
		 * Multply the Subtrahend by a scalar, without modifying it, then subtract from this vector.
		 */
		void multSub(double scalar, const snlVector& Subtrahend);

		void print();

		// The components of the vector.
		double components[4];
};

// Inline Functions
// ----------------

inline void snlVector::assign(const snlVector& copyFrom)
{
	components[0] = copyFrom.components[0];
	components[1] = copyFrom.components[1];
	components[2] = copyFrom.components[2];
	components[3] = copyFrom.components[3];
}

inline void snlVector::diff(const snlVector& v1, const snlVector& v2)
{
	components[0] = v2.components[0] - v1.components[0];
	components[1] = v2.components[1] - v1.components[1];
	components[2] = v2.components[2] - v1.components[2];
	components[3] = v2.components[3] - v1.components[3];
}

inline void snlVector::crossProduct (snlVector& v1, snlVector& v2)
{
	components[0] = (v1.components[1] * v2.components[2]) - (v1.components[2] * v2.components[1]);
	components[1] = (v1.components[2] * v2.components[0]) - (v1.components[0] * v2.components[2]);
	components[2] = (v1.components[0] * v2.components[1]) - (v1.components[1] * v2.components[0]);
	components[3] = 0.0;
}

inline double snlVector::dot (snlVector& vect)
{
	return (components[0] * vect.components[0] +
			components[1] * vect.components[1] +
			components[2] * vect.components[2] +
			components[3] * vect.components[3]);
}

inline double snlVector::lengthSqrd()
{
	return components[0] * components[0] +
		   components[1] * components[1] +
		   components[2] * components[2] +
		   components[3] * components[3];
}

inline void snlVector::projectXZ()
{
	components[1] = 0.0;  // y = 0.
}

inline void snlVector::projectXY()
{
	components[2] = 0.0;  // z = 0.
}

inline void snlVector::projectYZ()
{
	components[0] = 0.0;
}

inline double snlVector::x()
{
	return components[0];
}

inline double snlVector::y()
{
	return components[1];
}

inline double snlVector::z()
{
	return components[2];
}

inline double snlVector::w()
{
	return components[3];
}

inline void snlVector::x(double val)
{
	components[0] = val;
}

inline void snlVector::y(double val)
{
	components[1] = val;
}

inline void snlVector::z(double val)
{
	components[2] = val;
}

inline void snlVector::w(double val)
{
	components[3] = val;
}

inline void snlVector::setComponents(double x, double y, double z, double w)
{
	components[0] = x;
	components[1] = y;
	components[2] = z;
	components[3] = w;
}

inline void snlVector::setComponents(double x, double y, double z)
{
	components[0] = x;
	components[1] = y;
	components[2] = z;
	components[3] = 0.0;
}

inline void snlVector::getComponents(double* x, double* y, double* z, double* w)
{
	*x = components[0];
	*y = components[1];
	*z = components[2];
	*w = components[3];
}

inline void snlVector::getComponents(double* x, double* y, double* z)
{
	*x = components[0];
	*y = components[1];
	*z = components[2];
}

inline void snlVector::zero()
{
	components[0] = 0.0;
	components[1] = 0.0;
	components[2] = 0.0;
	components[3] = 0.0;
}

inline bool snlVector::isZero()
{

	return !(
		components[0] ||
		components[1] ||
		components[2] ||
		components[3]
	);
}

inline snlVector snlVector::operator * (double scalar) const
{
	return snlVector(
		components[0] * scalar,
		components[1] * scalar,
		components[2] * scalar,
		components[3] * scalar);
}

inline snlVector snlVector::operator * (const snlVector& vect) const
{
	// Duplicated from crossProduct function for speed.
	// [this] cross [vect].

	return snlVector(
		(components[1] * vect.components[2]) - (components[2] * vect.components[1]),
		(components[2] * vect.components[0]) - (components[0] * vect.components[2]),
		(components[0] * vect.components[1]) - (components[1] * vect.components[0]),
		0.0);
}

inline snlVector snlVector::operator + (const snlVector& vect) const
{
	return	snlVector(
		components[0] + vect.components[0],
		components[1] + vect.components[1],
		components[2] + vect.components[2],
		components[3] + vect.components[3]);
}

inline snlVector snlVector::operator - (const snlVector& vect) const
{
	return	snlVector(
		components[0] - vect.components[0],
		components[1] - vect.components[1],
		components[2] - vect.components[2],
		components[3] - vect.components[3]);
}

inline snlVector& snlVector::operator += (const snlVector& vect)
{
	components[0] += vect.components[0];
	components[1] += vect.components[1];
	components[2] += vect.components[2];
	components[3] += vect.components[3];

	return *this;
}

inline snlVector&  snlVector::operator -= (const snlVector& vect)
{
	components[0] -= vect.components[0];
	components[1] -= vect.components[1];
	components[2] -= vect.components[2];
	components[3] -= vect.components[3];

	return *this;
}

inline snlVector&  snlVector::operator *= (double scalar)
{
	components[0] *= scalar;
	components[1] *= scalar;
	components[2] *= scalar;
	components[3] *= scalar;

	return *this;
}

inline snlVector&  snlVector::operator /= (double scalar)
{
	components[0] /= scalar;
	components[1] /= scalar;
	components[2] /= scalar;
	components[3] /= scalar;

	return *this;
}

inline bool snlVector::operator == (snlVector& compare)
{
	return
		components[0] == compare.components[0] &&
		components[1] == compare.components[1] &&
		components[2] == compare.components[2] &&
		components[3] == compare.components[3];
}

inline void snlVector::multAdd(double scalar, const snlVector& addend)
{
	components[0] += scalar * addend.components[0];
	components[1] += scalar * addend.components[1];
	components[2] += scalar * addend.components[2];
	components[3] += scalar * addend.components[3];
}

inline void snlVector::multSub(double scalar, const snlVector& Subtrahend)
{
	components[0] -= scalar * Subtrahend.components[0];
	components[1] -= scalar * Subtrahend.components[1];
	components[2] -= scalar * Subtrahend.components[2];
	components[3] -= scalar * Subtrahend.components[3];
}

#endif
