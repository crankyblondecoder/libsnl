#ifndef SNLCTRLPOINTNET_H
#define SNLCTRLPOINTNET_H

#include "snlCtrlPoint.h"
#include "snlTransform.h"

#include <iostream>
#include <cmath>
using namespace std;

/**
 * Control Point Network - Base Class.
 */
class snlCtrlPointNet
{
	// Abstract base class of all control point objects.

	public:

		snlCtrlPointNet();
		virtual ~snlCtrlPointNet();

		snlCtrlPointNet(const snlCtrlPointNet&);  // Copy constructor.

		virtual void transform(unsigned ptIndex, snlTransform& transform);  // Transform a control point.

		virtual void transform(snlTransform& transform);  // Transform all points.

		virtual void transformSelected(snlTransform& transform);

		// Get number of control points that object holds.
		virtual unsigned getNumPts() const;

		// Get pointer to array of control points.
		virtual const snlCtrlPoint* getCtrlPts() const;

		// Get pointer to array of control points. Non constant return.
		virtual snlCtrlPoint* getCtrlPtsPtr();

		// Return copy of control point.
		snlCtrlPoint getPoint(unsigned index);

		// Return pointer to control point.
		virtual const snlCtrlPoint* getPointPtr(unsigned index);

		// Get transformed Z. Does not modify control points.
		virtual double getTransfZ(unsigned index, snlTransform&);

		virtual double getMaxTransfZ(snlTransform&);

		virtual double getMinTransfZ(snlTransform&);

		virtual bool hasPointsSelected();

		virtual unsigned numPointsSelected();

		virtual unsigned* getSelectedIndexes();

		virtual void selectAllPoints(bool yesNo = true);

		void selectPoint(unsigned index, bool yesNo = true);

		virtual bool isSelected(unsigned index);

		virtual void clearSelected();

		static double calcFlatness(snlPoint** points, unsigned size);
		double calcDeg1Flatness(snlPoint** points) const;

		double calcCurvature(snlPoint** points);

		bool isConvex(snlPoint** points, int numPts, double sensitivity = 0.0);

		virtual void replacePoints(snlCtrlPoint* newPoints);  // Replace control points with new ones.
		virtual void replacePoints(const snlCtrlPoint* newPoints, unsigned numNewPoints, unsigned replaceIndex,
			unsigned numToReplace);

		virtual void appendPointSpace(unsigned numPoints);  // Add space to end of array.
		virtual void truncatePointSpace(unsigned numPoints);  // Remove space from end of array.

		virtual void appendPoints(const snlCtrlPoint* points, unsigned numPoints);

		void print();
		void print(unsigned fromIndex, unsigned toIndex);

		bool hasConcurrentPoints() const;

		// Operators

		void operator +=(const snlCtrlPointNet& ctrlPointNet);
		void operator -=(const snlCtrlPointNet& ctrlPointNet);

		// *** Abstract Interface ***

		// Return list of connected control points.
		virtual int getCnctPts(unsigned index, snlCtrlPoint* retPts) = 0;

		// Return maximum number of connections a single control point can have.
		virtual int maxConnections() const = 0;

	protected:

		virtual bool checkBounds(unsigned index);  // Check index against array bounds.

		snlCtrlPoint* ctrlPts;
		unsigned ctrlPtArraySize;
};

#endif
