#ifndef SNL_CURVE_H
#define SNL_CURVE_H

#include "snlCtrlPointNetCurve.h"
#include "snlCurveBase.h"
#include "snlKnotVector.h"
#include "snlPoint.h"
#include "snlVector.h"
#include "snlVertex.h"
#include "snlVertexNet.h"

/**
 * General NURBS Curve
 */
class snlCurve : public snlCurveBase
{
	public:

		virtual ~snlCurve();

		snlCurve();

		/**
		 * Generate new curve that is a straight line between two points.
		 * @param deg Degree of curve.
		 * @param size Number of control points curve has.
		 * @param start Starting point of curve.
		 * @param end Ending point of curve.
		 */
		snlCurve(int degree, unsigned size, snlPoint& start, snlPoint& end);

		/**
		 * Generate new curve.
		 * @note Does not copy points and knots. So don't delete them elsewhere.
		 * @param deg Degree of curve.
		 * @param size Number of control points curve has.
		 * @param points Points to use.
		 * @param knots Knots to use.
		 */
		snlCurve(int degree, unsigned size, snlCtrlPoint* points, knot* knots = 0);

		/**
		 * Generate new curve.
		 * @note Does not copy points and knot vector. So don't delete them elsewhere.
		 * @param size Size of control point array.
		 * @param points Control point array.
		 * @param knotVect Knot vector to use.
		 */
		snlCurve(unsigned size, snlCtrlPoint* points, snlKnotVector* knotVector);

		/** Copy constructor. */
		snlCurve(const snlCurve& copyFrom);

		snlCurve& operator=(const snlCurve& curveToCopy);

		/**
		 * Interpolated / approximated curve.
		 * @note Array returned via retParams should be deleted by calling function.
		 * @param points Points to interpolate between.
		 * @param size Number of points.
		 * @param fittingType Type of interpolation from SNL_FITTING_TYPES.
		 * @param degree Resultant curve should be this degree.
		 * @param closedLoop The points specify a closed loop that should join smoothly.
		 * @param retParams Pointer to pointer that points to array of parameters that correspond to given points.
		 */
		snlCurve(snlPoint* points, unsigned size, int fittingType, int degree, bool closedLoop = false, knot** retParams = 0);

		/**
		 * Create curve as circular arc.
		 * @note Can not do an arc bigger than 180 degrees.
		 * @param startPoint Starting point of arc.
		 * @param endPoint Ending point of arc.
		 * @param centrePoint Centre point of arc.
		 * @param numSections Specify number of sections arc has. Use -1 to indicate an internally generated reasonable default
		 *        should be used.
		 */
		snlCurve(snlPoint& startPoint, snlPoint& endPoint, snlPoint& centrePoint, int numSections = - 1);

		/** Return reference to control point network object for curve. */
		snlCtrlPointNetCurve& controlPointNet();

		/**
		 * Evaluate Non Rational Homogeneous Curve Point.
		 * @param param Parameter to evaluate.
		 * @returns Homogeneous point on curve.
		 */
		snlPoint evalHmg(knot param) const;

		/**
		 * Evaluate rational non-homogeneous curve point.
		 * @param param Parameter to evaluate.
		 * @returns Non-homogeneous point on curve.
		 */
		virtual snlPoint eval(knot param) const;

		/**
		 * Evaluate Non Rational Homogeneous Surface Derivatives.
		 * @param param Parameter to evaluate at.
		 * @param deriv Derivative order to evaluate.
		 * @returns Array of snlPoint[deriv + 1]. Calling function must delete[] this array.
		*/
		snlPoint* evalDerivsHmg(knot param, unsigned deriv) const;

		/**
		 * Evaluate Rational Non-Homogeneous Surface Derivatives
		 * @param param Parameter to evaluate at.
		 * @param deriv Derivative order to evaluate.
		 * @returns Array of snlPoint[deriv + 1]. Calling function must delete[] array.
		 */
		snlPoint* evalDerivs(knot param, unsigned deriv) const;

		/**
		 * Velocity(first derivative) of curve.
		 * @param param Parameter to get velocity at.
		 */
		snlVector velocity(knot param);

		/**
		 * Insert a knot into knot vector and calculate new control points.
		 * @note ctrlPts MUST have an additional point space allocated at the end of each line in the array for the chosen
		 *       direction.
		 * @param iParam Parameter value to insert.
		 * @param reallocate Reallocate memory for control points.
		 */
		void insertKnot(knot iParam, bool reallocate);

		/**
		 * Insert multiple knots.
		 * @param iParam Parameter to insert.
		 * @param numToInsert Number of knots to insert.
		 * @param reallocate Reallocate memory for control points.
		 */
		void insertKnots(knot iParam, int numToInsert, bool reallocate);

		/**
		 * Remove multiple knots from index.
		 * @note Only removes multiples of the same parameter value initially at removal index.
		 * @param numKnots Number of knots to remove.
		 * @param removalIndex Index to remove knot from.
		 * @param tolerance Maximum error allowed before knot removal aborted. No tolerance if equals 0.
		 * @returns Tolerance achieved during knot removal whether successful or not.
		 */
		double removeKnots(int numKnots, unsigned removalIndex, double tolerance);

		/**
		 * Remove knot from curve.
		 * @param removalIndex Index to remove knot from.
		 * @param tolerance Maximum error allowed before knot removal aborted. No tolerance if equals 0.
		 * @returns Tolerance achieved during knot removal whether successful or not.
		 */
		double removeKnot(unsigned removalIndex, double tolerance);

		/** Refine control point net until tolerance is achieved. */
		void refine(double tolerance);

		/** Return maximum parameter value for curve. */
		double maxParam() const;

		/** Return minimum parameter value for curve. */
		double minParam() const;

		/** Return parameter at specified knot index. */
		double param(unsigned index) const;

		const snlKnotVector& knotVector();

		int degree();

		int size();

		/**
		 * Truncate curve.
		 * @param param Parameter to truncate at.
		 * @param keepLast Keep last part of curve instead of first part.
		 * @param reparam Reparameterise curve to original pararemeter boundaries.
		 */
		void truncate(knot param, bool keepLast = false, bool reparameterise = false);

		/**
		 * Insert partition into curve.
		 * @note Function basically makes sure degree knots are present at supplied parameter.
		 * @param param Parameter to insert partition into.
		 */
		void insertPartition(knot param);

		/**
		 * Do a linear Reparameterise on curve.
		 * @note Linear reparameterisations don't effect control points.
		 * @param startKnot New starting knot value of knot vector.
		 * @param endKnot Ending knot value of knot vector.
		 */
		void reparameterise(knot startKnot, knot endKnot);

		/** Reverse curves parametric evaluation direction. */
		void reverseEvalDirection();

		/**
		 * Global interpolation as closed loop.
		 * @note Array returned via retParams should be deleted by calling function.
		 * @param type Type of global interpolation from SNL_FITTING_TYPES.
		 * @param points Points to interpolate between.
		 * @param size Number of points.
		 * @param degree Resultant curve should be this degree.
		 * @param retParams Pointer to pointer that points to array of parameters that correspond to given points.
		 */
		void globalInterpClosedLoop(int type, snlPoint* points, unsigned size, int degree, knot** retParams);

		/**
		 * Global interpolation.
		 * @note Array returned via retParams should be deleted by calling function.
		 * @param type Type of global interpolation from SNL_FITTING_TYPES.
		 * @param points Points to interpolate between.
		 * @param size Number of points.
		 * @param degree Resultant curve should be this degree.
		 * @param retParams Pointer to pointer that points to array of parameters that correspond to given points.
		 */
		void globalInterp(int type, snlPoint* points, unsigned size, int degree, knot** retParams);

		/**
		 * Generate control points for global interpolation.
		 * @param points Array of points to interpolate between.
		 * @param size Size of points array.
		 * @param params Parameters that correspond to points array.
		 * @param knots Knot vector to use.
		 * @returns Array of control points the same size as the given "points" array.
		 */
		static snlCtrlPoint* genGlobalInterpPoints(snlPoint* points, unsigned size, knot* params, snlKnotVector* knots);

		/** Local quadratic interpolation. */
		void localInterpQuadratic(snlPoint* points, unsigned size);  // Local quadratic interpolation.

		/** Local cubic interpolation. */
		void localInterpCubic(snlPoint* points, unsigned size);  // Local cubic interpolation.

		/**
		 * Synchronise this curves knot vector to given curve.
		 * @param curve Curve to snync to.
		 * @note Knots are only ever added NOT removed. So if curve has less multiplicity at a particular span index then true
		 *       synchronisation will not occur and the caller should call the synchronise function on curve with this object
		 *       as it's argument.
		 */
		void synchronise(snlCurve& curve);

		/**
		 * Make this curve and given curve compatible.
		 * @param curve Curve to make this curve compatible with.
		 */
		void makeCompatible(snlCurve* curve);


		/**
		 * Create Bezier Segments over entire curve.
		 * @param retNumKnotsAdded Pointer to pointer to return array with number of inserted knots in it. Caller must delete
		 *        this array. First index in array corresponds to second knot span.
		 * @returns Number of elements in returned array.
		 */
		unsigned createBezierSegments(int** retNumKnotsAdded);

		/**
		 * Calculate new control points for Bezier segment that is being degree elevated.
		 * @param origDegree Original degree of segment.
		 * @param byDegree Number of degrees segment is being elevated by.
		 * @param origPoints Original points of non-elevated segment. Expected size is origDegree + 1.
		 * @param newPoints Array to store new points in. Must be size origDegree + byDegree + 1.
		 */
		static void elevateBezierSegmentPointsDegree(int origDegree, int byDegree, const snlCtrlPoint* origPoints,
			snlCtrlPoint* newPoints);

		/**
		 * Elevate degree of curve
		 * @param byDegree Number of degrees to elevate by.
		 */
		void elevateDegree(int byDegree);

		/**
		 * Append a curve to this curve.
		 * @param curveToAppend Curve to append to this curve.
		 * @param copy Copy curveToAppend and don't modify it in any way.
		 * @note This curve may have it's degree elevated to accomodate new curve.
		 */
		void appendCurve(snlCurve* curve, bool copy = true);

		/** Print curve to standard out. */
		void print();

		// Impl
		virtual void vertexNet(snlVertexNet* vNet, double tolerance, bool parametric);

		// Enumerations.

		enum SNL_FITTING_TYPES
		{
			SNL_GLOBAL_INTERPOLATION,
			SNL_GLOBAL_INTERP_CENTRIFUGAL,
			SNL_LOCAL_INTERPOLATION
		};

	protected:

		/**
		 * Generate vertex net based on evaluation of curve.
		 * @param vNet Vertex network to load with points.
		 * @param tolerance Tolerance to actual surface.
		 */
		void _vertexNetParam(snlVertexNet* vNet, double tolerance);

	private:

		int _deg; // Degree of curve.

		snlCtrlPointNetCurve* _ctrlPtNet;

		snlKnotVector* _knotVect;
};

#endif
