#ifndef SNL_CTRLPTNETCURVE_H
#define SNL_CTRLPTNETCURVE_H

#include "snlCtrlPointNet.h"

class snlCtrlPointNetCurve : public snlCtrlPointNet
{
	public:

		snlCtrlPointNetCurve(snlCtrlPoint* cPtArray, unsigned size, bool copy = false);

		snlCtrlPointNetCurve(unsigned size, snlPoint& start, snlPoint& end);

		virtual ~snlCtrlPointNetCurve();

		unsigned size() const;

		// Increase the control point net's size.
		snlCtrlPoint* grow();

		// Decrease the control point net's size.
		snlCtrlPoint* shrink();

		double calcFlatness(int index, int numPoints) const;

		void truncate(int atIndex, bool keepLast);

		void reverse();

		// snlCtrlPointNet Abstract implementation.

		virtual int getCnctPts(unsigned index, snlCtrlPoint* retPts);

		virtual int maxConnections() const;

	private:

};

#endif
