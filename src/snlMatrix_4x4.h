#ifndef SNL_MATRIX_4X4_H
#define SNL_MATRIX_4X4_H

#include "snlVector.h"

#include <iostream>

using namespace std;

/** 4x4 Matrix of doubles */
class snlMatrix_4X4
{
	public:

		/** Constructs a new matrix that is set to identity. */
		snlMatrix_4X4();

		/** Constructs a new matrix that is initialised to given elements. */
		snlMatrix_4X4(double initialElements[16]);

		/** Set matrix to identity. */
		void ident();

		/** Setup matrix as translation identity. */
		void translateIdent(double x, double y, double z);

		/**
		 * Create rotation identity for rotation about x axis.
		 * @param yy y coordinate scalar for new y.
		 * @param yz z coordinate scalar for new y.
		 * @param zy y coordinate scalar for new z.
		 * @param zz z coordinate scalar for new z.
		 */
		void rotateXIdent(double yy, double yz, double zy, double zz);

		/**
		 * Create rotation identity for rotation about y axis.
		 * @param xx x coordinate scalar for new x.
		 * @param xz z coordinate scalar for new x.
		 * @param zx x coordinate scalar for new z.
		 * @param zz z coordinate scalar for new z.
		 */
		void rotateYIdent(double xx, double xz, double zx, double zz);

		/**
		 * Create rotation identity for rotation about z axis.
		 * @param xx x coordinate scalar for new x.
		 * @param xy y coordinate scalar for new x.
		 * @param yx x coordinate scalar for new y.
		 * @param yy y coordinate scalar for new y.
		 */
		void rotateZIdent(double xx, double xy, double yx, double yy);

		/**
		 * Create scaling identity for scaling operation.
		 * @param x Scaling factor in X direction.
		 * @param y Scaling factor in Y direction.
		 * @param z Scaling factor in Z direction.
		 */
		void scaleIdent(double x, double y, double z);

		/**
		 * Pre-multiply the given matrix to this matrix and store result in this.
		 * @note Pre-multiplication covers all standard use cases for 3D scene graphs. The higher transforms always accumulate
		 *       by pre-multiplying as the graph is traversed.
		 * @param multMatrix Matrix to multiply to this.
		 */
		void preMultiply(snlMatrix_4X4&);

		/**
		 * Multiply this matrix to one vector. Storing the result in the vector.
		 * This is the matrix multiplied by the vector treated as a column vector.
		 */
		void multiply(snlVector& vector);

		/**
		 * Multiply this matrix to one or more vectors. Storing the result in each vector.
		 * This is the matrix multiplied by the vector treated as a column vector.
		 * @param vectors Pointer to first element in array of vectors to process.
		 * @param numVectors Number of vectors of array to process.
		 */
		void multiply(snlVector* vectors, unsigned numVectors);

		void print();  //!< Print matrice to standard out.

	protected:

		/** Matrix elements in Column Major order (ie OpenGL Standard). */
		double _elements[16];
};

inline void snlMatrix_4X4::ident()
{
	_elements[0] = 1.0;
	_elements[1] = 0.0;
	_elements[2] = 0.0;
	_elements[3] = 0.0;

	_elements[4] = 0.0;
	_elements[5] = 1.0;
	_elements[6] = 0.0;
	_elements[7] = 0.0;

	_elements[8] = 0.0;
	_elements[9] = 0.0;
	_elements[10] = 1.0;
	_elements[11] = 0.0;

	_elements[12] = 0.0;
	_elements[13] = 0.0;
	_elements[14] = 0.0;
	_elements[15] = 1.0;
}

inline void snlMatrix_4X4::translateIdent(double x, double y, double z)
{
	_elements[0] = 1.0;
	_elements[1] = 0.0;
	_elements[2] = 0.0;
	_elements[3] = 0.0;

	_elements[4] = 0.0;
	_elements[5] = 1.0;
	_elements[6] = 0.0;
	_elements[7] = 0.0;

	_elements[8] = 0.0;
	_elements[9] = 0.0;
	_elements[10] = 1.0;
	_elements[11] = 0.0;

	_elements[12] = x;
	_elements[13] = y;
	_elements[14] = z;
	_elements[15] = 1.0;
}

inline void snlMatrix_4X4::rotateXIdent(double yy, double yz, double zy, double zz)
{
	_elements[0] = 1.0;
	_elements[1] = 0.0;
	_elements[2] = 0.0;
	_elements[3] = 0.0;

	_elements[4] = 0.0;
	_elements[5] = yy;
	_elements[6] = zy;
	_elements[7] = 0.0;

	_elements[8] = 0.0;
	_elements[9] = yz;
	_elements[10] = zz;
	_elements[11] = 0.0;

	_elements[12] = 0.0;
	_elements[13] = 0.0;
	_elements[14] = 0.0;
	_elements[15] = 1.0;
}

inline void snlMatrix_4X4::rotateYIdent(double xx, double xz, double zx, double zz)
{
	_elements[0] = xx;
	_elements[1] = 0.0;
	_elements[2] = zx;
	_elements[3] = 0.0;

	_elements[4] = 0.0;
	_elements[5] = 1.0;
	_elements[6] = 0.0;
	_elements[7] = 0.0;

	_elements[8] = xz;
	_elements[9] = 0.0;
	_elements[10] = zz;
	_elements[11] = 0.0;

	_elements[12] = 0.0;
	_elements[13] = 0.0;
	_elements[14] = 0.0;
	_elements[15] = 1.0;
}

inline void snlMatrix_4X4::rotateZIdent(double xx, double xy, double yx, double yy)
{
	_elements[0] = xx;
	_elements[1] = yx;
	_elements[2] = 0.0;
	_elements[3] = 0.0;

	_elements[4] = xy;
	_elements[5] = yy;
	_elements[6] = 0.0;
	_elements[7] = 0.0;

	_elements[8] = 0.0;
	_elements[9] = 0.0;
	_elements[10] = 1.0;
	_elements[11] = 0.0;

	_elements[12] = 0.0;
	_elements[13] = 0.0;
	_elements[14] = 0.0;
	_elements[15] = 1.0;
}

inline void snlMatrix_4X4::scaleIdent(double x, double y, double z)
{
	_elements[0] = x;
	_elements[1] = 0.0;
	_elements[2] = 0.0;
	_elements[3] = 0.0;

	_elements[4] = 0.0;
	_elements[5] = y;
	_elements[6] = 0.0;
	_elements[7] = 0.0;

	_elements[8] = 0.0;
	_elements[9] = 0.0;
	_elements[10] = z;
	_elements[11] = 0.0;

	_elements[12] = 0.0;
	_elements[13] = 0.0;
	_elements[14] = 0.0;
	_elements[15] = 1.0;
}

#endif
