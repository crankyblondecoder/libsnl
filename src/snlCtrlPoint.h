#ifndef SNLCTRLPOINT_H
#define SNLCTRLPOINT_H

#include "snlPoint.h"

class snlCtrlPoint : public snlPoint
{
	public:

		snlCtrlPoint();

		snlCtrlPoint(snlPoint&);

		void operator =(const snlPoint& copyFrom);

		void select(bool yesNo);
		bool isSelected();

		void weight(double setTo);
		double weight();

	private:

		bool selected; // Point has been selected.
};

#endif
