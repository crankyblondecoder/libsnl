// *** Curve Base Class ***

#ifndef SNL_CURVE_BASE_H
#define SNL_CURVE_BASE_H

#include "snlKnotVector.h"
#include "snlPoint.h"
#include "snlVertexNet.h"

class snlCurveBase
{
	public:

		virtual ~snlCurveBase(){};

		virtual snlPoint eval(knot param) const = 0;

		/**
		 * Return approximation to curve.
		 * @param tolerance Tolerance to approximate to.
		 * @param parametric Do a parametric analysis as opposed to knot refinement.
		 * @param vNet Vertex net to fill with data.
		 */
		virtual void vertexNet(snlVertexNet* vNet, double tolerance, bool parametric) = 0;
};

#endif
