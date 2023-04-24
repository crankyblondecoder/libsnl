// *** Curve Base Class ***

#ifndef SNL_CURVE_BASE_H
#define SNL_CURVE_BASE_H

#include "snlKnotVector.h"
#include "snlPoint.h"
#include "snlVertexNet.h"

#ifdef SGI_MIPS

    #include <iostream.h>
    #include <math.h>
    #include <float.h>

#else

    #include <iostream>
    #include <cmath>
    #include <cfloat>

    using namespace std;

#endif

class snlCurveBase
{
    public:

        virtual ~snlCurveBase(){};

        virtual snlPoint eval ( knot param ) const = 0;

        virtual void vertexNet ( snlVertexNet* vNet, double tolerance, bool parametric ) = 0;
};

#endif
