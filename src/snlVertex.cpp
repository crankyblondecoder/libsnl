#include "snlVertex.h"

snlVertex::snlVertex()
{
}

snlVertex::snlVertex(double x, double y, double z, double w)
	: snlPoint(x, y, z, w)
{
}

void snlVertex::normal(snlVector& setTo)
{
	_norm = setTo;
}

snlVector& snlVertex::normal()
{
	return _norm;
}

void snlVertex::evalParamU(knot value)
{
	_paramU = value;
}

knot snlVertex::evalParamU()
{
	return _paramU;
}

void snlVertex::evalParamV(knot value)
{
	_paramV = value;
}

knot snlVertex::evalParamV()
{
	return _paramV;
}
