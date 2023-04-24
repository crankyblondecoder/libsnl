// *** Base Class Of All Surfaces ***

#ifndef SNL_SURFACE_BASE_H
#define SNL_SURFACE_BASE_H

#include "snlKnotVector.h"
#include "snlPoint.h"
#include "snlTriangleMesh.h"
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

class snlSurfaceBase
{
    public:

        virtual ~snlSurfaceBase(){};

        virtual snlPoint eval ( knot paramU, knot paramV ) const = 0;

        virtual void vertexNet ( snlVertexNet* vNet, double tolerance, bool parametric ) = 0;

        virtual void triangleMesh ( snlTriangleMesh* triMesh, int toleranceType, double tolerance ) = 0;

        enum meshToleranceType
        {
            SNL_TOL_DISTANCE,  // Must be within distance tolerance to surface.
            SNL_TOL_ANGLE      // Angle between successive sections must be less than tolerance.
        };
};

#endif
