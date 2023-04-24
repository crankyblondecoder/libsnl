// *** Surface of Revolution ***

#include "snlSurfaceOfRevolution.h"
#include "snlUtil.h"

snlSurfaceOfRevolution::snlSurfaceOfRevolution()
{
    profile = 0;
    axis_start = 0;
    axis_end = 0;

    rot_angle = M_PI * 2.0;
}

snlSurfaceOfRevolution::~snlSurfaceOfRevolution()
{
    if ( profile ) delete profile;
    if ( axis_start ) delete axis_start;
    if ( axis_end ) delete axis_end;
}

snlSurfaceOfRevolution::snlSurfaceOfRevolution ( snlCurve* profileCurve, snlPoint* axisStart, snlPoint* axisEnd, double rotationAngle )
{
    // Construct new surface of revolution.
    // ------------------------------------
    // profileCurve:    Curve to revolve about axis.
    // axisStart:       Start of revolution axis.
    // axisEnd:         End of revolution axis.
    // rotationAngle:   Angle to revolve about axis.
    //
    // Notes:    All objects are owned by this object.
    //           Revolution about axis is right handed.

    profile = profileCurve;
    axis_start = axisStart;
    axis_end = axisEnd;

    if ( rotationAngle > M_PI * 2.0 || rotationAngle == 0.0)
        rot_angle = M_PI * 2.0;
    else if ( rotationAngle < ( - M_PI * 2.0 ) )
        rot_angle = - M_PI * 2.0;
    else
        rot_angle = rotationAngle;
}

snlCurve& snlSurfaceOfRevolution::profileCurve()
{
    return * profile;
}

void snlSurfaceOfRevolution::profileCurve ( snlCurve* profileCurve )
{
    if ( profile ) delete profile;

    profile = profileCurve;
}

snlPoint& snlSurfaceOfRevolution::axisStart()
{
    return *axis_start;
}

void snlSurfaceOfRevolution::axisStart ( snlPoint* startPoint )
{
    if ( axis_start ) delete axis_start;

    axis_start = startPoint;
}

snlPoint& snlSurfaceOfRevolution::axisEnd()
{
    return *axis_end;
}

void snlSurfaceOfRevolution::axisEnd ( snlPoint* endPoint )
{
    if ( axis_end ) delete axis_end;

    axis_end = endPoint;
}

double snlSurfaceOfRevolution::rotationAngle()
{
    return rot_angle;
}

void snlSurfaceOfRevolution::rotationAngle ( double angle )
{
    rot_angle = angle;
}

void snlSurfaceOfRevolution::vertexNet ( snlVertexNet* vNet, double tolerance, bool parametric )
{
    // Return approximation to surface.
    // --------------------------------
    // tolerance:    Tolerance to approximate to.
    // parametric:   Do a parametric analysis as opposed to knot refinement.
    // vNet:         Vertex net to fill with data.

    snlCurve* profileCopy = new snlCurve ( *profile );

    if ( tolerance > 0.0 )
        profileCopy -> refine ( tolerance );

    const snlCtrlPoint* ctrlPts = profileCopy -> controlPointNet().getCtrlPts();

    int numPts = profileCopy -> controlPointNet().size();

    vNet -> vertexNet ( ctrlPts, numPts );

    // Find angle step to use.

    snlPoint axis_start_norm ( *axis_start );
    axis_start_norm.project();

    snlPoint axis_end_norm ( *axis_end );
    axis_end_norm.project();

    // Get largest radius

    double maxRadius = 0.0;

    for ( int index = 0; index < numPts; index ++ )
    {
        double dist = distToLine ( axis_start_norm, axis_end_norm, ctrlPts [ index ] );

        if ( dist > maxRadius )
            maxRadius = dist;
    }

    // Calculate steps based on tolerance and maximum radius.

    double angleStep = 2 * acos ( 1.0 - ( tolerance / maxRadius ) );

    int numSteps = (int ) ( rot_angle / angleStep ) + 1;

    angleStep = rot_angle / (double) numSteps;

    // Rotate and append points at discrete angle steps.

    snlTransform transf;

    transf.rotate ( angleStep, axis_start_norm, axis_end_norm );

    for ( int step = 0; step < numSteps; step ++ )
    {
        profileCopy -> controlPointNet().transform ( transf );

        vNet -> appendRow ( ctrlPts );
    }
}

snlPoint snlSurfaceOfRevolution::eval ( knot paramU, knot paramV ) const
{
    snlPoint retPoint;

    // !@#$ Incomplete

    return retPoint;
}

void snlSurfaceOfRevolution::triangleMesh ( snlTriangleMesh* triMesh, double tolerance )
{
}

