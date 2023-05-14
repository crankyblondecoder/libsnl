#include "snlCtrlPoint.h"

snlCtrlPoint::snlCtrlPoint()
{
	selected = false;
}

snlCtrlPoint::snlCtrlPoint(snlPoint& pt)
	: snlPoint(pt)
{
	selected = false;
}

void snlCtrlPoint::operator = (const snlPoint& copyFrom)
{
	// Copy data from another point.
	// -----------------------------

	components[0] = copyFrom.components[0];
	components[1] = copyFrom.components[1];
	components[2] = copyFrom.components[2];
	components[3] = copyFrom.components[3];

	selected = false;
}

void snlCtrlPoint::select(bool yesNo)
{
	// Set selection state of control point.
	// -------------------------------------

	selected = yesNo;
}

bool snlCtrlPoint::isSelected()
{
	return selected;
}

void snlCtrlPoint::weight(double setTo)
{
	// Set control points weight.
	// --------------------------

	if(components[3] == 0.0)
	{
		components[3] = setTo;
	}
	else
	{
		double multFactor = setTo / components[3];

		for(int index = 0; index < 4; index ++)
			components[index] *= multFactor;
	}
}

double snlCtrlPoint::weight()
{
	return components[3];
}

