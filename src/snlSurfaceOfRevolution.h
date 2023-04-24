// ***!! Deprecated - Use snlSurface Instead !!****

// *** Surface of Revolution ***

#ifndef SNL_SURFACEOFREVOLUTION_H
#define SNL_SURFACEOFREVOLUTION_H

#include "snlCurve.h"
#include "snlPoint.h"
#include "snlSurfaceBase.h"

class snlSurfaceOfRevolution : public snlSurfaceBase
{
    public:

        snlSurfaceOfRevolution();
        virtual ~snlSurfaceOfRevolution();

        snlSurfaceOfRevolution ( snlCurve* profileCurve, snlPoint* axisStart, snlPoint* axisEnd, double rotationAngle );

        snlCurve& profileCurve();
        void profileCurve ( snlCurve* profileCurve );

        snlPoint& axisStart();
        void axisStart ( snlPoint* startPoint );

        snlPoint& axisEnd();
        void axisEnd ( snlPoint* endPoint );

        double rotationAngle();
        void rotationAngle ( double angle );

        // Abstract Implementation.

        virtual void vertexNet ( snlVertexNet* vNet, double tolerance, bool parametric );

        virtual snlPoint eval ( knot paramU, knot paramV ) const;

        virtual void triangleMesh ( snlTriangleMesh* triMesh, double tolerance );

    private:

        snlCurve*    profile;

        snlPoint*    axis_start;
        snlPoint*    axis_end;

        double       rot_angle;
};

#endif
