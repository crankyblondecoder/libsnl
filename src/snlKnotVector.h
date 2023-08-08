#ifndef SNLKNOTVECTOR_H
#define SNLKNOTVECTOR_H

#include <iostream>
#include <cmath>

using namespace std;

// Limiting the size of the degree allows memory optimisations.
#define SNL_KNOT_VECTOR_MAX_DEG 7
#define SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS 8

typedef double knot;
typedef double basis;

/**
 * Single dimension array of knots.
 */
class snlKnotVector
{
	public:

		virtual ~snlKnotVector();

		/// Copy constructor.
		snlKnotVector(const snlKnotVector&);

		/**
		 * Construct using an existing knot array.
		 * @param knotArrayToUse Array of knots.
		 * @param size Size of array
		 * @param degree Degree of knot vector.
		 * @param knotVectType Type of knot vector.
		 * @param copy If true copy the knot array.
		 */
		snlKnotVector(knot* knotArrayToUse, unsigned size, int degree, int knotVectType = 1, bool copy = false);

		/**
		 * Generate a new Knot Vector
		 * @note Assumes clamped (open) vector.
		 * @param startVal Starting parametric knot value.
		 * @param endVal End parametric knot value.
		 * @param numKnots Number of knots in vector.
		 * @param degree Degree to use during evaluation of basis functions.
		*/
		snlKnotVector(knot startVal, knot endVal, unsigned numKnots, int degree);

		/**
		 * Generate knot vector given existing parameters.
		 * @note Used for interpolation.
		 * @param params Parameters to average between.
		 * @param size Size of parameters array.
		 * @param degree Degree of vector.
		*/
		snlKnotVector(int size, int degree, knot* params);

		/**
		 * Generate knot vector for Bezier patch.
		 * @param degree Degree of knot vector.
		 */
		snlKnotVector(int degree);

		/**
		 * Assignment Operator.
		 */
		snlKnotVector& operator=(const snlKnotVector& KnotVectToCopy);

		enum knotVectorType
		{
			/// Clamped
			open = 1,
			/// Periodic will not be supported for some time yet :-)
			periodic = 2
		};

		knot getKnotVal(unsigned index) const;

		const knot* getKnotPtr(unsigned index);

		/// Get the size of knot vector.
		unsigned getSize() const;

		/// Get degree of knot vector.
		int getDegree();

		/// Set degree of knot vector.
		void setDegree(int degree);

		bool equals(const snlKnotVector& knotVect) const;

		/** Insert new knot into vector.
		 * @param param knot to insert.
		 * @param numTimes Number of times to insert knot into vector.
		 */
		void insertKnot(knot param, int numTimes = 1);

		/// Remove knot at spanIndex
		void removeKnot(unsigned spanIndex);

		/**
		 * Increase size of knot vector array.
		 * @param bySize Size to increase knot vector by.
		 */
		void grow(unsigned bySize);

		/** Increase multiplicity of knots at span.
		 * @param spanIndex Index of knot span to process.
		 * @param numKnotsToAdd Number of knots to add at spanIndex.
		 */
		void increaseMultiplicity(unsigned spanIndex, int numKnotsToAdd);

		/**
		 * Find Knot Span corresponding to parameter
		 * @param param Parameter to find span of.
		 */
		unsigned findSpan(knot param) const;

		/// Find the number of non-zero length spans.
		unsigned numSpans() const;

		// Find the knot index of first non-zero span.
		unsigned firstSpan() const;

		// Find next non-zero length span given spanIndex to start from.
		unsigned nextSpan(unsigned spanIndex) const;

		// Find previous non-zero length span given spanIndex to start from.
		unsigned previousSpan(unsigned spanIndex) const;

		/**
		 * Find the knot multiplicity at index
		 * @param index Index to search.
		 * @returns Number of knots found.
		 */
		int findMultiplicity(unsigned index) const;

		/**
		 * Find the knot multiplicity at a particular knot value.
		 * @param param Parameter to evaluate.
		 * @returns Returns Number of knots found.
		 */
		int findMultiplicity(knot param) const;

		/**
		 * Truncate knot vector.
		 * @note Truncates at last of knots valued at param.
		 * @note Assumes degree knots are present at param.
		 * @param param Parameter to truncate at.
		 * @param keepLast Keep last section of knot vector, discard first part.
		 */
		void truncate(knot param, bool keepLast);

		/**
		 * Reparameterise knot vector using a linear reparameterisation.
		 * @param startKnot Knot vectors new starting value.
		 * @param endKnot Knot vectors new ending value.
		 */
		void reparameterise(knot startKnot, knot endKnot);

		/// Reverse knot vector.
		void reverse();

		/**
		 * Join another knot vector to the end of this one.
		 * @note Assumes end clamp of this and start clamp values of other knot vector are the same. This is _not_ just a
		 *       join of arrays, it results in a properly formed knot vector.
		 * @param knotVector Knot vector to join. Must be of same degree as this knot vector.
		 */
		void join(snlKnotVector* knotVector);

		// Get pointer to the array of knots contained in this.
		const knot* getKnotArray();

		/// Find maximum knot value in vector.
		knot max() const;
		/// Find minimum knot value in vector.
		knot min() const;

		// Get the knot vector type.
		int getType() const;

		/**
		 * Evaluate Basis Functions
		 * @param param Paramater value to process at.
		 * @param retArray Array to return evaluated values in. Must be size deg + 1.
		 * @returns Array of basis function values, evaluated at param.
		 */
		void evalBasis(knot param, basis* retArray);

		/**
		 * Evaluate basis functions and their derivatives.
		 * @param param Parameter to process at.
		 * @param deriv Which derivative to process to.
		 * @param retArray Two dimensional array of basis and basis derivative values. Must be size [deriv + 1][deg + 1].
		 *        ie Each row of the array holds the evaluated basis function values for that particular derivative.
		 *        TODO describe which deriviative appears first in array ...
		 */
		void evalBasisDeriv(knot param, int deriv, basis* retArray);

		/**
		 * Calculate alphas used for knot removal.
		 * @param span Span of knot where removal is taking place.
		 * @returns Array of alphas. Must be deleted by caller.
		 */
		double* calcRemovalAlphas(unsigned span);

		/// Print knot vector to std out.
		void print();

		/**
		 * Print knot vector to std out.
		 * @note Prints in c++ code format.
		 */
		void print_cpp();

	private:

		void __copyFrom(const snlKnotVector& vector);

		knot* _knots;
		unsigned _vectorSize;

		int _deg; // Degree associated with vector.

		// Type of knot vector
		int _kvType;
};

#endif
