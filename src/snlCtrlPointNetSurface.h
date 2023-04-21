#ifndef SNLCTRLPTNETSURFACE_H
#define SNLCTRLPTNETSURFACE_H

#include "snlCtrlPointNet.h"

class snlCtrlPointNetSurface : public snlCtrlPointNet
{
	public:

		snlCtrlPointNetSurface(snlCtrlPoint* cPtArray, unsigned sizeU, unsigned sizeV, bool copy = false);

		snlCtrlPointNetSurface(unsigned sizeU, unsigned sizeV, snlPoint& origin,
			snlPoint& cornerU, snlPoint& cornerV);

		virtual ~snlCtrlPointNetSurface();

		unsigned getSizeU() const;
		unsigned getSizeV() const;

		void setSizeU(unsigned size);
		void setSizeV(unsigned size);

		snlCtrlPoint* getPoint(unsigned indexU, unsigned indexV) const;

		// Increase the control point net's size.
		snlCtrlPoint* growU(int increaseBy = 1, bool reallocate = true);
		snlCtrlPoint* growV(int increaseBy = 1, bool reallocate = true);

		// Decrease the control point net's size.
		snlCtrlPoint* shrinkU();
		snlCtrlPoint* shrinkV();

		double calcFlatness(int indexU, int indexV, int numPointsU, int numPointsV);

		double calcFlatnessU(int indexU, int indexV, int numPoints, bool degree1) const;
		double calcFlatnessV(int indexU, int indexV, int numPoints, bool degree1) const;

		double maxFlatnessU(int span);
		double maxFlatnessV(int span);

		double maxCurvatureU();
		double maxCurvatureV();

		void locatePoints(int indexU, int indexV, int numPointsU, int numPointsV, snlCtrlPoint** pointsLocated) const;

		void locatePointsU(int indexU, int indexV, int numPoints, snlCtrlPoint** testPoints) const;
		void locatePointsV(int indexU, int indexV, int numPoints, snlCtrlPoint** testPoints) const;

		void selectPoint(int indexU, int indexV);
		void selectLineConstU(int indexU); // Select line in constant U direction.
		void selectLineConstV(int indexV); // Select line in constant V direction.

		void print();
		void printCompare(snlCtrlPointNetSurface& compareTo);
		void print_cpp();

		// snlCtrlPointNet Abstract implementation.

		virtual int getCnctPts(unsigned index, snlCtrlPoint* retPts);

		virtual int maxConnections() const;

	private:

		int sizeU, sizeV;  // Array size in U and V directions.

};

#endif


