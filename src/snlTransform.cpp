#include "snlTransform.h"

#include <cmath>
using namespace std;

snlTransform::snlTransform()
	: snlMatrix_4X4()
{
}

snlTransform::snlTransform(double initialElements[16])
	: snlMatrix_4X4(initialElements)
{
}

void snlTransform::translate(double x, double y, double z)
{
	if(x == 0.0 && y == 0.0 && z == 0.0) return;

	__scratchMatrix.translateIdent(x, y, z);

	preMultiply(__scratchMatrix);
}

void snlTransform::rotateX(double angle)
{
	__scratchMatrix.rotateXIdent(cos(angle), - sin(angle), sin(angle), cos(angle));

	preMultiply(__scratchMatrix);
}

void snlTransform::rotateY(double angle)
{
	__scratchMatrix.rotateYIdent(cos(angle), sin(angle), - sin(angle), cos(angle));

	preMultiply(__scratchMatrix);
}

void snlTransform::rotateZ(double angle)
{
	__scratchMatrix.rotateZIdent(cos(angle), - sin(angle), sin(angle), cos(angle));

	preMultiply(__scratchMatrix);
}

void snlTransform::rotate(double angle, snlPoint& axisStart, snlPoint& axisEnd)
{
	snlVector axisDir = axisEnd - axisStart;

	rotate(angle, axisStart, axisDir);
}

void snlTransform::rotate(double angle, snlPoint& axisStart, snlVector& axisDirection)
{
	axisDirection.normalise();

	// Move axis to origin
	translate(- axisStart.x(), - axisStart.y(), - axisStart.z());

	if(axisDirection.y() == 0.0 && axisDirection.z() == 0.0)
	{
		// Just a rotation about the x axis.

		if(axisDirection.x() > 0.0)
		{
			rotateX(angle);
		}
		else
		{
			rotateX(- angle);
		}
	}
	else if(axisDirection.x() == 0.0 && axisDirection.z() == 0.0)
	{
		// Just a rotation about the Y axis.

		if(axisDirection.y() > 0.0)
		{
			rotateY(angle);
		}
		else
		{
			rotateY(- angle);
		}
	}
	else if(axisDirection.x() == 0.0 && axisDirection.y() == 0.0)
	{
		// Just a rotation about the x axis.

		if(axisDirection.z() > 0.0)
		{
			rotateZ(angle);
		}
		else
		{
			rotateZ(- angle);
		}
	}
	else
	{
		// Rotate about x axis to yz plane.
		double projHyp = sqrt(axisDirection.y() * axisDirection.y() + axisDirection.z() * axisDirection.z());

		double sinRotX = axisDirection.y() / projHyp;
		double cosRotX = axisDirection.z() / projHyp;

		__scratchMatrix.rotateXIdent(cosRotX, - sinRotX, sinRotX, cosRotX);
		preMultiply(__scratchMatrix);

		// Rotate about y axis to z axis.
		double sinRotY = -(axisDirection.x());
		double cosRotY = projHyp;

		__scratchMatrix.rotateYIdent(cosRotY, sinRotY, - sinRotY, cosRotY);
		preMultiply(__scratchMatrix);

		// Rotate about z axis.
		rotateZ(angle);

		// Inverse rotation about y axis. ie Rotate back into yz plane. Negate angle of rotation.
		__scratchMatrix.rotateYIdent(cosRotY, - sinRotY, sinRotY, cosRotY);
		preMultiply(__scratchMatrix);

		// Inverse rotation about x axis.
		__scratchMatrix.rotateXIdent(cosRotX, sinRotX, - sinRotX, cosRotX);
		preMultiply(__scratchMatrix);
	}

	// Inverse translation.
	translate(axisStart.x(), axisStart.y(), axisStart.z());
}

void snlTransform::scale(double x, double y, double z)
{
	__scratchMatrix.scaleIdent(x, y , z);

	preMultiply(__scratchMatrix);
}

void snlTransform::align(snlVector& vector1, snlVector& vector2)
{
	// Get angle to rotate through.
	double rotAngle = vector1.angle(vector2);

	if(rotAngle == 0.0) return;

	// Generate normal to both vectors.
	snlVector normal;
	normal.crossProduct(vector1, vector2);

	snlTransform tmpTransform;

	snlPoint origin;

	// Rotate about normal.
	if((rotAngle != 0.0) &&(!normal.isZero()))
	{
		tmpTransform.rotate(rotAngle, origin, normal);
	}

	preMultiply(tmpTransform);
}
