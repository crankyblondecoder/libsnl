#ifndef SNL_CIRCULAROFFSET_CURVE_H
#define SNL_CIRCULAROFFSET_CURVE_H

#include "snlCurve.h"
#include "snlCurveBase.h"
#include "snlPoint.h"

/**
 * NURBS Curve of Circular Offset from Base Curve
 */
class snlCircularOffsetCurve : public snlCurveBase
{
	public:

		snlCircularOffsetCurve();
		virtual ~snlCircularOffsetCurve();

		snlCircularOffsetCurve(snlCircularOffsetCurve& copyFrom);

		snlCircularOffsetCurve(snlCurve* baseCurve, snlPoint* axisStart, snlPoint* axisEnd);

		void refine(double tolerance);

		virtual void applyOffset(snlPoint& point, snlPoint chordOffset, snlPoint angleOffset, snlPoint tangentOffset) const;

		int numOffsets();
		void generateOffsets(int type, double startOffset, double endOffset);
		void offset(int index, int type, double val, double weight = 1.0);
		double offset(int index, int type);

		void vertexNetParam(snlVertexNet* vNet, double tolerance);

		int size();

		double maxParam() const;
		double minParam() const;

		enum offsetType
		{
			CHORD,
			ANGLE,
			TANGENT
		};

		// Abstract Implementation.

		virtual void vertexNet(snlVertexNet* vNet, double tolerance, bool parametric);

		virtual snlPoint eval(knot param) const;

	protected:

		snlCurve* base_curve;

		snlCurve* chord_offsetCurve;  // Offsets that correspond to control points.
		snlCurve* angle_offsetCurve;
		snlCurve* tangent_offsetCurve;

		snlPoint* axis_start;
		snlPoint* axis_end;
};

#endif
