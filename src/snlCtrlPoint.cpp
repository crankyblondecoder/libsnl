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

	elements[0] = copyFrom.elements[0];
	elements[1] = copyFrom.elements[1];
	elements[2] = copyFrom.elements[2];
	elements[3] = copyFrom.elements[3];

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

	if(elements[3] == 0.0)
	{
		elements[3] = setTo;
	}
	else
	{
		double multFactor = setTo / elements[3];

		for(int index = 0; index < 4; index ++)
			elements[index] *= multFactor;
	}
}

double snlCtrlPoint::weight()
{
	return elements[3];
}

