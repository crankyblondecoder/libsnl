#ifndef SNLTRANSFORM_H
#define SNLTRANSFORM_H

#include "snlPoint.h"
#include "snlVector.h"
#include "snlMatrix_4x4.h"

/**
 * Transformation defined in R⁴.
 * @note For efficiency of calculation transforms only support pre-multiplication.
 * @note To transform a vector use the underlying superclass multiply function.
 */
class snlTransform : public snlMatrix_4X4
{
	public:

		snlTransform();

		snlTransform(double initialElements[16]);

		/**
		 * Apply translation.
		 * @param x Translation in x direction.
		 * @param y Translation in y direction.
		 * @param z Translation in z direction.
		 */
		void translate(double x, double y, double z);

		/**
		 * Rotate about +ve x axis
		 * @param angle Rotation angle in radians.
		 */
		void rotateX(double angle);

		/**
		 * Rotate about +ve y axis
		 * @param angle Rotation angle in radians.
		 */
		void rotateY(double angle);

		/**
		 * Rotate about +ve z axis
		 * @param angle Rotation angle in radians.
		 */
		void rotateZ(double angle);

		/**
		 * Rotate about an arbitrary axis.
		 * @param angle Angle to rotate about axis in radians.
		 * @param axisStart Start of axis.
	     * @param axisDirection Direction of axis. Warning: This vector will be normalised.
		 */
		void rotate(double angle, snlPoint& axisStart, snlVector& axisDirection);

		/**
		 * Rotate about an arbitrary axis.
		 * @param angle Angle to rotate about axis in radians.
		 * @param axisStart Start point of line that describes the axis.
	     * @param axisEnd End point of line that describes the axis.
		 */
		void rotate(double angle, snlPoint& axisStart, snlPoint& axisEnd);

		/** Apply scaling factors. */
		void scale(double x, double y, double z);

		/**
		 * Generate transform that will align vector1 to vector2.
		 * Effectively rotates vector1 about normal to plane described by vector1 and vector2.
		 * @param vector1 Vector to align.
		 * @param vector2 Vector to align vector 1 to.
		 * @note Can be used for transforming coordinate systems.
		*/
		void align(snlVector& vector1, snlVector& vector2);

	private:

		/** Used for adding operations to the transform. */
		snlMatrix_4X4 _scratchMatrix;
};

#endif
