// *** Test program for libSNL ***

#include <iostream>
#include <cmath>

using namespace std;

#include "../snlSquareLinear.h"
#include "../snlCurve.h"
#include "../snlSurface.h"
#include "../snlCtrlPoint.h"
#include "../snlVector.h"

snlCurve* generateCurve()
{
   // Create some points to interpolate.

	snlPoint* points = new snlPoint[10];

	double x = 0.0;
	double y = 0.0;
	double z = 10.0;
	double w = 1.0;

	double zIncr = 1.0;

	for(int count = 0; count < 10; count ++)
	{
		points[count].components(x, y, z, w);

		x += 1.0;
		y += 1.0;
		z -= zIncr;

		zIncr -= 0.2;
	}

	knot* params;

	// Generate curve.

	snlCurve* curve = new snlCurve(points, 10, snlCurve::SNL_GLOBAL_INTERPOLATION, 3, false, &params);

	delete[] params;
	delete[] points;

	return curve;
}

snlCurve* generateReverseCurve()
{
   // Create some points to interpolate.

	snlPoint* points = new snlPoint[10];

	double x = 0.0;
	double y = 0.0;
	double z = 10.0;
	double w = 1.0;

	double zIncr = 1.0;

	for(int count = 9; count >= 0; count --)
	{
		points[count].components(x, y, z, w);

		x += 1.0;
		y += 1.0;
		z -= zIncr;

		zIncr -= 0.2;
	}

	knot* params;

	// Generate curve.

	snlCurve* curve = new snlCurve(points, 10, snlCurve::SNL_GLOBAL_INTERPOLATION, 3, false, &params);

	delete[] points;
	delete[] params;

	return curve;
}

snlSurface* generateSurface()
{
	// Generate surface for use elsewhere in the test program.
	// -------------------------------------------------------

	snlCtrlPoint* points = new snlCtrlPoint[100];

	double x = - 5.0;
	double y = 0.0;
	double z = 10.0;
	double w = 1.0;

	double xyGrad = 0.02;
	double xyGradIncr = 0.01;

	double xIncr = 0.02;
	double yIncr = 0.05;
	double zIncr = -0.2;

	double xStep = 1.0;

	for(int tIndex = 0; tIndex < 10; tIndex ++)
	{
		for(int uIndex = 0; uIndex < 10; uIndex ++)
		{
			double newX = x +(xStep * uIndex);
			double newY = y +(newX - x) * xyGrad;

			points[( tIndex * 10) + uIndex].components(newX, newY, z, w);

			x += xIncr;
			y += yIncr;
			z += zIncr;

			xyGrad += xyGradIncr;
		}
	}

	snlSurface* retSurface = new snlSurface(5, 5, 10, 10, points, 0, 0);

	return retSurface;
}

snlSurface* generateSurface2()
{
	// Generate surface for use elsewhere in the test program.
	// -------------------------------------------------------

	snlCtrlPoint* points = new snlCtrlPoint[110];

	double z = 2.5;
	double w = 1.0;

	points[0].components(- 2.5, 0.5, z, w);
	points[1].components(- 2.0, 1.0 , z, w);
	points[2].components(- 1.5, 0.5, z, w);
	points[3].components(- 1.0, 0.0, z, w);
	points[4].components(- 0.5, 0.5, z, w);
	points[5].components(0.0, 1.0, z, w);
	points[6].components(0.5, 0.5, z, w);
	points[7].components(1.0, 0.0, z, w);
	points[8].components(1.5, 0.5, z, w);
	points[9].components(2.0, 1.0, z, w);
	points[10].components(2.5, 1.5, z, w);

	double zIncr = - 0.5;

	double angleStep = 0.314;

	snlTransform rotation;
	snlTransform translation;

	snlPoint axisStart(0.0, 0.0, 0.0);
	snlPoint axisEnd(0.0, 0.0, 1.0);

	rotation.rotate(angleStep, axisStart, axisEnd);
	translation.translate(0.0, 0.0, zIncr);

	int index = 11;

	for(int uIndex = 1; uIndex < 10; uIndex ++)
	{
		for(int vIndex = 0; vIndex < 11; vIndex ++)
		{
			points[index] = points[vIndex];

			rotation.multiply(points[index]);
			translation.multiply(points[index]);

			index ++;
		}

		rotation.rotate(angleStep, axisStart, axisEnd);
		translation.translate(0.0, 0.0, zIncr);
	}

	snlSurface* retSurface = new snlSurface(5, 5, 10, 11, points, 0, 0);

	return retSurface;
}

snlSurface* generateSawToothSurface()
{
	// Generate surface for use elsewhere in the test program.
	// -------------------------------------------------------
	// Notes:   Generates a surface suitable for concave polygon detection.

	int size = 11;
	int deg = 5;

	snlCtrlPoint* points = new snlCtrlPoint[size * size];

	double x = - 5.0;
	double y = 0.0;
	double z = 5.0;
	double w = 1.0;

	double ySaw = 0.5;

	double xStep = 1.0;
	double yIncr = 0.1;
	double zIncr = -1.0;

	for(int tIndex = 0; tIndex < size; tIndex ++)
	{
		for(int uIndex = 0; uIndex < size; uIndex ++)
		{
			points[( tIndex * size) + uIndex].components(x +(xStep * uIndex), y +(uIndex % 2) * ySaw, z, w);
		}

		y += yIncr;
		z += zIncr;
	}

	knot* knotsV = new knot[size + deg + 1];

	for(int count = 0; count < deg + 1; count ++) knotsV[count] = 0.0;
	for(int count = deg + 1; count < deg * 2 + 1; count ++) knotsV[count] = 0.5;
	for(int count = deg * 2 + 1; count < size + deg + 1; count ++) knotsV[count] = 1.0;

	snlSurface* retSurface = new snlSurface(5, 5, 11, 11, points, 0, knotsV);

	return retSurface;
}

snlSurface* generateSurfaceOfRevolution()
{
	// Generate surface for use elsewhere in the test program.
	// -------------------------------------------------------

	// Code used from vTurbine.

	double impeller_depth = 2.0;
	double inlet_radius = 1.0;
	double outlet_radius = 3.0;
	double innerProfile_weight1 = 0.5;
	double innerProfile_weight2 = 0.5;

	snlPoint* start = new snlPoint(0.0, inlet_radius, impeller_depth, 1.0);
	snlPoint* end = new snlPoint(0.0, outlet_radius, 0.0, 1.0);  // Axis is rooted at origin.

	snlCurve innerProfileCurve(3, 6, *start, *end);

	snlCtrlPoint* cpts = innerProfileCurve.controlPointNet().getCtrlPtsPtr();

	cpts[1].components(0.0, inlet_radius,(impeller_depth * 2.0) / 3.0, 1.0);
	cpts[2].components(0.0, inlet_radius, impeller_depth / 3.0, 1.0);
	cpts[3].components(0.0,(outlet_radius - inlet_radius) / 3.0 + inlet_radius, 0.0, 1.0);
	cpts[4].components(0.0,(( outlet_radius - inlet_radius) * 2.0) / 3.0 + inlet_radius, 0.0, 1.0);


	cpts[2].weight(innerProfile_weight1);
	cpts[3].weight(innerProfile_weight2);

	start -> components(0.0, 0.0, impeller_depth, 1.0);
	end -> components(0.0, 0.0, 0.0, 1.0);

	return new snlSurface(innerProfileCurve, *start, *end, 270.0);
}

snlSurface* generateAmbigSurface()
{
	// Generate surface for ambiguous edge detection.

	int degree = 3;

	snlCtrlPoint* line = new snlCtrlPoint[degree + 1];

	line[0].components(0.0, 0.0, 5.0);
	line[1].components(0.0, 1.0, 1.0);
	line[2].components(0.0, - 1.0, 1.0);
	line[3].components(0.0, 0.0, 5.0);

	snlCtrlPoint* ctrlPts = new snlCtrlPoint[( degree + 1) *(degree + 1)];

	double xAdjust = 1.0;

	int index = 0;

	for(int uIndex = 0; uIndex < degree + 1; uIndex ++)
	{
		for(int vIndex = 0; vIndex < degree + 1; vIndex ++)
		{
			ctrlPts[index] = line[vIndex];
			ctrlPts[index].x(ctrlPts[index].x() + xAdjust * uIndex);

			index ++;
		}
	}

	delete[] line;

	return new snlSurface(degree, degree, degree + 1, degree + 1, ctrlPts, 0, 0);
}

snlSurface* generateAmbigSurface2()
{
	// Generate surface for ambiguous edge detection.

	int degree = 3;

	snlCtrlPoint* line = new snlCtrlPoint[degree + 1];

	line[0].components(0.0, 0.0, 5.0);
	line[1].components(0.0, 1.0, 1.0);
	line[2].components(0.0, - 1.0, 1.0);
	line[3].components(0.0, 0.0, 5.0);

	snlCtrlPoint* ctrlPts = new snlCtrlPoint[( degree + 1) *(degree + 1)];

	double xAdjust = 1.0;

	int index = 0;

	for(int uIndex = 0; uIndex < degree + 1; uIndex ++)
	{
		for(int vIndex = 0; vIndex < degree + 1; vIndex ++)
		{
			ctrlPts[index] = line[uIndex];
			ctrlPts[index].x(ctrlPts[index].x() + xAdjust * vIndex);

			index ++;
		}
	}

	delete[] line;

	return new snlSurface(degree, degree, degree + 1, degree + 1, ctrlPts, 0, 0);
}

bool test_surfaceInterp()
{
	bool passed = true;

	double tolerance = 0.0001;

	// Generate surface from existing surfaces control points.

	snlSurface* testSurf = generateSurface();

	const snlCtrlPoint* ctrlPts = testSurf -> controlPoints();

	int numU = testSurf -> sizeU();
	int numV = testSurf -> sizeV();

	int numPts = numU * numV;

	snlPoint* pointsInterp = new snlPoint[numPts];

	for(int index = 0; index < numPts; index ++)
		pointsInterp[index] = ctrlPts[index];

	snlSurface* interpSurf = new snlSurface(snlSurface::SNL_GLOBAL_INTERP_CENTRIFUGAL, pointsInterp,
											  numU, numV, 4, 4);

	// Project orginal data points to surface and calculate error.

	double maxError = 0.0;

	snlPoint* toProject = new snlPoint[numPts];

	for(int index = 0; index < numPts; index ++)
		toProject[index] = pointsInterp[index];

	int numProjLocns;

	snlSurfLocn* projLocns = interpSurf -> fastProject(toProject, numPts, &numProjLocns, tolerance, 0.00001, 2, 1, 1);

	for(int index = 0; index < numPts; index ++)
	{
		if(maxError < projLocns[index].dist) maxError = projLocns[index].dist;
	}

	delete[] projLocns;
	delete[] toProject;

	if(maxError > tolerance) passed = false;

	cout << "Max Error: " << maxError << " - ";

	delete testSurf;
	delete interpSurf;
	delete[] pointsInterp;

	if(! passed)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	return passed;
}

bool test_surfaceRefine()
{
	bool passed = true;

	double tolerance = 0.05;

	snlSurface* testSurf = generateSurface();
	snlSurface* testSurf2 = generateSurface2();

	testSurf -> refine(tolerance);
	//testSurf2 -> refineHullBezier(tolerance);  // Broken with no time to fix.
	testSurf2 -> refine(tolerance);

	// Project all control points to surface and calculate error.

	int numSurfaces = 1;

	snlSurface* surfs[2];

	surfs[0] = testSurf;
	surfs[1] = testSurf2;

	double maxError = 0.0;

	cout << "\n\n";

	cout << "\tFinished refinement, starting projections. This may take a while.\n\n";

	for(int surfNum = 0; surfNum < numSurfaces; surfNum ++)
	{
		snlSurface* cSurf = surfs[surfNum];

		int numPoints =(cSurf -> controlPointNet()).getNumPts();

		cout << "\tTesting Surface Number: " << surfNum << "  Number of projections: " << numPoints << "\n";

		const snlCtrlPoint* ctrlPoints = cSurf -> controlPoints();

		snlPoint* toProject = new snlPoint[numPoints];

		for(int index = 0; index < numPoints; index ++)
			toProject[index] = ctrlPoints[index];

		int numProjLocns;

		snlSurfLocn* projLocns = cSurf -> fastProject(toProject, numPoints, &numProjLocns, tolerance, 0.00001, 2, 1, 1);

		for(int index = 0; index < numPoints; index ++)
		{
			if(maxError < projLocns[index].dist) maxError = projLocns[index].dist;
		}

		delete[] projLocns;
		delete[] toProject;
	}

	if(maxError > tolerance) passed = false;

	cout << "\tMax Error: " << maxError << " - ";

	delete testSurf;
	delete testSurf2;

	if(! passed)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	return passed;
}

bool test_ambig()
{
	bool passed = true;

	snlSurface* testSurf = generateAmbigSurface();
	snlSurface* testSurf2 = generateAmbigSurface2();

	sEdge surfEdges[4];
	sEdge surfEdges2[4];

	int numEdges = testSurf -> hasAmbigEdges(surfEdges, 1.0e-6);
	int numEdges2 = testSurf2 -> hasAmbigEdges(surfEdges2, 1.0e-6);

	if(numEdges == 0) passed = false;
	if(numEdges2 == 0) passed = false;

	for(int index = 0; index < numEdges; index ++)
	{
		if(surfEdges[index].direction != 0) passed = false;
		if(surfEdges2[index].direction != 1) passed = false;
	}

	delete testSurf;
	delete testSurf2;

	if(! passed)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	return passed;
}

bool testSurfaceOfRevolution()
{
	// Test accuracy of surface of rotation.
	// -------------------------------------

	cout << "\n\n";

	bool passed = true;

	snlSurface* origSurf = generateSurfaceOfRevolution();
	snlSurface* testSurf = generateSurfaceOfRevolution();

	// Test knot insertion at multiple locations.

	int numSteps = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = minParamU + paramStepU;
	double paramV = minParamV;

	double maxError = 0;

	int index = 0;

	testSurf -> refineHull_U(0.005);
	testSurf -> refineHull_V(0.005);

	for(int indexU = 0; indexU < numSteps; indexU ++)
	{
		paramV = minParamV;

		for(int indexV = 0; indexV < numSteps; indexV ++)
		{

			paramV += paramStepV;

			if(paramV > maxParamV) paramV = maxParamV;

			index ++;
		}

		paramU += paramStepU;

		if(paramU > maxParamU) paramU = maxParamU;
	}

	// Clean up, report and return.

	cout << "Max Error = " << maxError << " - ";

	delete testSurf;
	delete origSurf;

	if(! passed)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	return passed;
}

bool testSurfacePointInversion()
{
	// Test Surface Point Inversion Function
	// -------------------------------------

	//cout << "\n\n";

	bool passed = true;

	snlSurface* testSurf = generateSurface();

	int numSteps = 10;

	int maxPass = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = 0.0;
	double paramV = 0.0;

	double maxError = 0;

	double iterTol = 1.0e-8;
	double normTol = 1.0e-6;

	int numEval = numSteps * numSteps;

	// Fill array with evaluated points.

	snlPoint* evalPts = new snlPoint[numEval];

	int index = 0;

	for(int indexU = 0; indexU < numSteps; indexU ++)
	{
		paramV = minParamV;

		for(int indexV = 0; indexV < numSteps; indexV ++)
		{
			evalPts[index] = testSurf -> eval(paramU, paramV);

			//cout << "(" << indexU << ", " << indexV << ") ";
			//cout << "EvalParamU: " << paramU << "  EvalParamV: " << paramV << "\n";

			paramV += paramStepV;

			if(paramV > maxParamV) paramV = maxParamV;

			index ++;
		}

		paramU += paramStepU;

		if(paramU > maxParamU) paramU = maxParamU;
	}

	int numReturned;

	snlSurfLocn* inverted = testSurf -> invert(evalPts, numEval, &numReturned, iterTol, normTol, maxPass);

	//cout << "Num Sent: " << numEval << " Num Returned: " << numReturned << "\n";

	for(index = 0; index < numReturned; index ++)
	{
		snlVector delta = inverted[index].pt - evalPts[inverted[index].origPtIndex];

		if(maxError < delta.length())
			maxError = delta.length();
	}

	delete[] inverted;
	delete[] evalPts;

	if(maxError > iterTol) passed = false;

	// Clean up, report and return.

	cout << "Tolerance = " << iterTol << ", Max Error = " << maxError << " - ";

	if(numEval > numReturned)
	{
		cout << "Number Sent: " << numEval << " Number Returned: " << numReturned << " - ";
		passed = false;
	}

	delete testSurf;

	if(! passed)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	return passed;
}

bool testSurfaceConvexBezierSegmentation()
{
	// Test decomposition of surface into Bezier segements.
	// ----------------------------------------------------

	cout << "\n\n";

	bool passed = true;

	snlSurface* origSurf = generateSawToothSurface();
	snlSurface* testSurf = generateSawToothSurface();

	// Test convex detection function.

	snlPoint testPoints[] = { snlPoint(- 3.0, 0.0, 0.0),
							  snlPoint(- 1.5, 2.5, 0.0),
							  snlPoint(1.5, - 2.0, 0.0),
							  snlPoint(2.5, 1.5, 0.0),
							  snlPoint(3.0, 2.5, 0.0)
							};

	snlPoint testPoints2[] = { snlPoint(- 3.0, 0.0, 0.0),
							   snlPoint(- 1.5, 2.5, 0.0),
							   snlPoint(1.5, 5.0, 0.0),
							   snlPoint(2.5, 1.5, 0.0),
							   snlPoint(3.0, 2.5, 0.0)
							 };


	snlPoint* testPointPtrs[] = { testPoints, testPoints + 1, testPoints + 2, testPoints + 3, testPoints + 4 };
	snlPoint* testPointPtrs2[] = { testPoints2, testPoints2 + 1, testPoints2 + 2, testPoints2 + 3, testPoints2 + 4 };

	if(( testSurf -> controlPointNet()).isConvex(testPointPtrs, 5)
		 &&  !(testSurf -> controlPointNet()).isConvex(testPointPtrs2, 5))
	{
		passed = false;
		cout << "\tConvex Points Test - Failed\n";
	}
	else
		cout << "\tConvex Points Test - Passed\n";

	// Test knot insertion at multiple locations.

	int numSteps = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = minParamU + paramStepU;
	double paramV = minParamV;

	double maxError = 0;

	int numV, numU;

	testSurf -> createConvexBezierSegments(&numU, &numV);

	//cout << "Num U Segments: " << numU << "  Num V Segments: " << numV << "\n";

	for(int indexU = 1; indexU < numSteps - 1; indexU ++)
	{
		paramV = minParamV + paramStepV;

		for(int indexV = 1; indexV < numSteps - 1; indexV ++)
		{
			snlPoint original = origSurf -> eval(paramU, paramV);

			snlPoint test = testSurf -> eval(paramU, paramV);

			double error = (test - original).length();

			if(error > maxError) maxError = error;

			paramV += paramStepV;
		}

		paramU += paramStepU;
	}

	delete origSurf;
	delete testSurf;

	if(maxError > 2.0e-13) passed = false;

	cout << "\tSurface Decomposition Into Convex Bezier Segments - Max Error = " << maxError << "\n";

	if(! passed)
		cout << "\tFailed\n";
	else
		cout << "\tPassed\n";

	return passed;
}

bool testSurfaceBezierSegmentation()
{
	// Test decomposition of surface into Bezier segements.
	// ----------------------------------------------------

	cout << "\n\n";

	bool passed = true;

	snlSurface* origSurf = generateSurface();
	snlSurface* testSurf = generateSurface();

	// Test knot insertion at multiple locations.

	int numSteps = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = minParamU + paramStepU;
	double paramV = minParamV;

	double maxError = 0;

	testSurf -> createBezierSegments();

	for(int indexU = 1; indexU < numSteps - 1; indexU ++)
	{
		paramV = minParamV + paramStepV;

		for(int indexV = 1; indexV < numSteps - 1; indexV ++)
		{
			snlPoint original = origSurf -> eval(paramU, paramV);

			snlPoint test = testSurf -> eval(paramU, paramV);

			double error = (test - original).length();

			if(error > maxError) maxError = error;

			paramV += paramStepV;
		}

		paramU += paramStepU;
	}

	delete origSurf;
	delete testSurf;

	if(maxError > 2.0e-13) passed = false;

	cout << "\tSurface Decomposition Into Bezier Segments - Max Error = " << maxError << "\n";

	if(! passed)
		cout << "\tFailed\n";
	else
		cout << "\tPassed\n";

	return passed;
}

bool testSurfaceKnotInsert()
{
	// Test knot insertion functions.
	// ------------------------------

	cout << "\n\n";

	bool passed = true;

	snlSurface* origSurf = generateSurface();
	snlSurface* testSurf = generateSurface();
	snlSurface* testSurf2 = generateSurface();

	// Test knot insertion at multiple locations.

	int numSteps = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = minParamU + paramStepU;
	double paramV = minParamV;

	double maxError = 0;
	double maxError2 = 0;

	// Warning!! Don't insert knots at end clamps.

	for(int indexU = 1; indexU < numSteps - 1; indexU ++)
	{
		paramV = minParamV + paramStepV;

		int multi =(origSurf -> knotVectorU()).findMultiplicity(paramU);

		testSurf -> insertKnot(paramU, snlSurface::SNL_U_DIR);

		testSurf2 -> insertKnot(paramU, snlSurface::SNL_U_DIR, testSurf2 -> degreeU() - multi);
		//testSurf2 -> insertKnot(paramU, snlSurface::SNL_U_DIR, 4);

		//( testSurf2 -> controlPointNet()).printCompare(testSurf -> controlPointNet());

		for(int indexV = 1; indexV < numSteps - 1; indexV ++)
		{
			// Evaluate before knot insertion.

			snlPoint original = origSurf -> eval(paramU, paramV);

			// Insert knot and evaluate.

			if(indexU < 2)
			{
				multi =(origSurf -> knotVectorV()).findMultiplicity(paramV);

				testSurf -> insertKnot(paramV, snlSurface::SNL_V_DIR);

				testSurf2 -> insertKnot(paramV, snlSurface::SNL_V_DIR, testSurf2 -> degreeV() - multi);
				//testSurf2 -> insertKnot(paramV, snlSurface::SNL_V_DIR, 4);
			}

			snlPoint inserted = testSurf -> eval(paramU, paramV);
			snlPoint inserted2 = testSurf2 -> eval(paramU, paramV);

			double error = (inserted - original).length();
			double error2 =(inserted2 - original).length();

			//cout << "Param U, V: " << paramU << ", " << paramV << "Error: " << error << " Error2: " << error2 << "\n";

			if(error > maxError) maxError = error;
			if(error2 > maxError2) maxError2 = error2;

			paramV += paramStepV;
		}

		paramU += paramStepU;
	}

	delete origSurf;
	delete testSurf;
	delete testSurf2;

	if(maxError > 1.5e-13) passed = false;
	if(maxError2 > 1.5e-13) passed = false;

	cout << "\tSurface Single Knot Insertion - Max Error = " << maxError << "\n";
	cout << "\tSurface Multiple Knot Insertion - Max Error = " << maxError2 << "\n";

	if(! passed)
		cout << "\tFailed\n";
	else
		cout << "\tPassed\n";

	return passed;
}

bool testDerivEval()
{
	// Test derivative evaluation.
	// ---------------------------

	cout << "\n\n";

	bool passed = true;

	snlSurface* testSurf = generateSurface();

	int numSteps = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = minParamU;
	double paramV = minParamV;

	double deltaU_val =(maxParamU - minParamU) / 1000000.0;
	double deltaV_val =(maxParamV - minParamV) / 1000000.0;

	double maxError = 0;
	double maxMixedPartialError = 0;

	for(int indexU = 0; indexU < numSteps; indexU ++)
	{
		paramV = minParamV;

		for(int indexV = 0; indexV < numSteps; indexV ++)
		{
			// Get derivatives
			snlPoint* derivs = testSurf -> evalDerivs(paramU, paramV, 2, 2);  // Calculate 2nd derivs as well.

			double deltaU, deltaV;

			if(paramU < maxParamU)
				deltaU = deltaU_val;
			else
				deltaU = - deltaU_val;

			if(paramV < maxParamV)
				deltaV = deltaV_val;
			else
				deltaV = - deltaV_val;

			snlVector velocityU(derivs[3]);
			snlVector velocityV(derivs[1]);

			// Get points a little tiny delta from evaluated point.

			// Vo.t + 0.5 a.(t*t).

			snlPoint deltaPointU = derivs[0];

			deltaPointU += velocityU * deltaU;

			deltaPointU += derivs[6] *(deltaU * deltaU * 0.5);

			snlPoint deltaPointV = derivs[0];

			deltaPointV += velocityV * deltaV;

			deltaPointV += derivs[2] *(deltaV * deltaV * 0.5);

			// Calculate distance between actual and approximated points.

			snlVector deltaVectorU = deltaPointU - testSurf -> eval(paramU + deltaU, paramV);
			snlVector deltaVectorV =  deltaPointV - testSurf -> eval(paramU, paramV + deltaV);

			double distU = deltaVectorU.length();
			double distV = deltaVectorV.length();

			//cout << "ParamU: " << paramU << "  ParamV: " << paramV << "  distU: " << distU << " distV: " << distV << "\n";

			// Asses first mixed partial.

			snlVector deltaVelocityU(derivs[4]);
			deltaVelocityU *= deltaV;
			deltaVelocityU += velocityU;

			snlVector deltaVelocityV(derivs[4]);
			deltaVelocityV *= deltaU;
			deltaVelocityV += velocityV;


			snlPoint* derivsDeltaU = testSurf -> evalDerivs(paramU + deltaU, paramV, 2, 2);
			snlPoint* derivsDeltaV = testSurf -> evalDerivs(paramU, paramV + deltaV, 2, 2);

			double discrMixedU = deltaVelocityU.dot(derivsDeltaV[3]) - snlVector(derivsDeltaV[3]).lengthSqrd();
			discrMixedU /= snlVector(derivsDeltaV[3]).lengthSqrd();

			double discrMixedV = deltaVelocityV.dot(derivsDeltaU[1]) - snlVector(derivsDeltaU[1]).lengthSqrd();
			discrMixedV /= snlVector(derivsDeltaU[1]).lengthSqrd();

			delete[] derivsDeltaU;
			delete[] derivsDeltaV;

			//cout << "ParamU: " << paramU << "  ParamV: " << paramV << "  discrU: " << discrMixedU << " discrV: " << discrMixedV << "\n";

			// Final processing

			paramV += paramStepV;

			if(paramV > maxParamV) paramV = maxParamV;

			if( distU > maxError) maxError = distU;

			if( distV > maxError) maxError = distV;

			if(discrMixedV > maxMixedPartialError) maxMixedPartialError = discrMixedV;
			if(discrMixedU > maxMixedPartialError) maxMixedPartialError = discrMixedU;

			delete[] derivs;
		}

		paramU += paramStepU;

		if(paramU > maxParamU) paramU = maxParamU;
	}

	if(maxError > 1.0e-12 || maxMixedPartialError > 1.0e-9) passed = false;

	// Clean up, report and return.

	cout << "\tSurface Derivatives, velocity / acceleration approximation - Max Error = " << maxError << "\n";
	cout << "\tSurface Derivatives, mixed partial velocity approximation - Max Error = " << maxMixedPartialError << "\n";

	delete testSurf;

	if(! passed)
	{
		cout << "\t\tFailed\n";
		return false;
	}
	else
		cout << "\tPassed\n";

	return passed;
}

bool testSurfaceProjection()
{
	// Test Surface Projection
	// -----------------------

	cout << "\n\n";

	bool passed = true;

	//snlSurface* testSurf = generateSurface();
	snlSurface* testSurf = generateSurface2();

	int numSteps = 10;

	int maxPass = 10;

	// Step through surface and evaluate.

	double minParamU =(testSurf -> knotVectorU()).min();
	double minParamV =(testSurf -> knotVectorV()).min();
	double maxParamU =(testSurf -> knotVectorU()).max();
	double maxParamV =(testSurf -> knotVectorV()).max();

	double paramStepU =(maxParamU - minParamU) /(numSteps - 1);
	double paramStepV =(maxParamV - minParamV) /(numSteps - 1);

	double paramU = 0.0;
	double paramV = 0.0;

	double maxError = 0;
	double maxErrorFast = 0;

	double normTol = 1.0e-6;
	double iterTol = 1.0e-8;

	double normLength = 0.01;

	int numEval = numSteps * numSteps;

	// Fill array with evaluated points.

	snlPoint* evalPts = new snlPoint[numEval];

	int index = 0;

	for(int indexU = 0; indexU < numSteps; indexU ++)
	{
		paramV = minParamV;

		for(int indexV = 0; indexV < numSteps; indexV ++)
		{
			snlPoint point;
			snlVector velU;
			snlVector velV;

			testSurf -> velocities(paramU, paramV, point, velU, velV);

			snlVector normal;

			normal.crossProduct(velU, velV);

			normal.length(normLength);

			point += normal;

			evalPts[index] = point;

/*
if(index == 1)
{
			//cout << "(" << indexU << ", " << indexV << ") ";
			cout << index << "- EvalParamU: " << paramU << "  EvalParamV: " << paramV << "\n";

	int numReturned;
	snlSurfLocn* projected = testSurf -> project(evalPts + index, 1, &numReturned, iterTol, normTol, maxPass);
	for(int retIndex = 0; retIndex < numReturned; retIndex ++)
	{
		cout << "Num Returned: " << numReturned << "\n";
		cout << "ParamU: " << projected[retIndex].paramU << " ParamV: " << projected[retIndex].paramV
			 << "  Dist: " << projected[retIndex].dist << "\n";
	}
	delete[] projected;
}
*/
			paramV += paramStepV;

			if(paramV > maxParamV) paramV = maxParamV;

			index ++;
		}

		paramU += paramStepU;

		if(paramU > maxParamU) paramU = maxParamU;
	}

	int numReturned;
	int numReturnedFast;

	snlSurfLocn* projected = testSurf -> project(evalPts, numEval, &numReturned, iterTol, normTol, maxPass);
	snlSurfLocn* fastProjected = testSurf -> fastProject(evalPts, numEval, &numReturnedFast, iterTol, normTol, maxPass,
														   1, 2);

	//cout << "Fast -> Num Sent: " << numEval << " Num Returned: " << numReturnedFast << "\n";

	for(index = 0; index < numReturned; index ++)
	{
		if(maxError < projected[index].cos)
			maxError = projected[index].cos;

		if(maxErrorFast < fastProjected[index].cos)
			maxErrorFast = fastProjected[index].cos;

		if(fastProjected[index].cos > normTol)
		{
			cout << "Fast - not under tolerance - ptIndex: " << fastProjected[index].origPtIndex << "\n";
		}

		if(index)
		{
			if(projected[index].origPtIndex - projected[index - 1].origPtIndex > 1)
				cout << "Hole Found - start: " <<  projected[index - 1].origPtIndex << "  end: "
					 << projected[index].origPtIndex << "\n";

			if(fastProjected[index].origPtIndex - fastProjected[index - 1].origPtIndex > 1)
				cout << "Hole Found - start: " <<  fastProjected[index - 1].origPtIndex << "  end: "
					 << fastProjected[index].origPtIndex << "\n";
		}
		else
		{
			if(projected[index].origPtIndex != 0)
				cout << "Starts at: " << projected[index].origPtIndex << "\n";

			if(fastProjected[index].origPtIndex != 0)
				cout << "Starts at: " << fastProjected[index].origPtIndex << "\n";
		}
	}

	delete[] projected;
	delete[] fastProjected;
	delete[] evalPts;

	if(maxError > normTol) passed = false;
	if(maxErrorFast > normTol) passed = false;

	// Clean up, report and return.

	cout << "\tTolerance = " << normTol << ", Max Error = " << maxError << " - "
		 << "Fast Max Error = " << maxErrorFast << " - ";

	if(numEval > numReturned)
	{
		cout << "\nNumber Sent: " << numEval << " Number Returned: " << numReturned << " - ";
		passed = false;
	}

	if(numEval > numReturnedFast)
	{
		cout << "\nFast Proj -> Number Sent: " << numEval << " Number Returned: " << numReturned << " - ";
		passed = false;
	}

	delete testSurf;

	if(! passed)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	return passed;
}

bool testSquareLinear()
{
   // Test square linear class.

	double* coeffs = new double[16];

	double coeffVals[16] = { -2,  2, -4,  -6,
								-3,  6,  3, -15,
								 5, -8, -1,  17,
								 1,  1, 11,   7 };

	for(int index = 0; index < 16; index ++)
		coeffs[index] = coeffVals[index];

	double * rhs = new double[12];

	double rhsVals[12] = { -4,  2, -7,
							  -3,  8,  4,
							   9, -2, -5,
							   7, 10,  6 };

	for(int index = 0; index < 12; index ++)
		rhs[index] = rhsVals[index];

	snlSquareLinear solver(4, 3, coeffs, rhs);

	//solver.print();

	// Compare the solutions to known values.

	float knownVals[12] = { 35, 655.0 / 6.0, - 1175.0 / 12.0,
							   28, 554.0 / 6.0, - 1006.0 / 12.0,
							   -7, -142.0 / 6.0, 278.0 / 12.0,
							   3, 59.0 / 6.0, -115.0 / 12.0 };

	bool success = true;

	for(int index = 0; index < 12; index ++)
	{
		if(( (float) rhs[index]) != knownVals[index]) success = false;
	}

	return success;
}

bool testCurveInterpolation()
{
	// Create some points to interpolate.

	snlPoint* points = new snlPoint[10];

	double x = 0.0;
	double y = 0.0;
	double z = 10.0;
	double w = 1.0;

	double zIncr = 1.0;

	for(int count = 0; count < 10; count ++)
	{
		points[count].components(x, y, z, w);

		x += 1.0;
		y += 1.0;
		z -= zIncr;

		zIncr -= 0.2;
	}

	knot* params;

	// Generate curve.

	cout << "\n\tGlobal Interpolation ... ";
	snlCurve curve(points, 10, snlCurve::SNL_GLOBAL_INTERPOLATION, 3, false, &params);

	double maxError = 0.0;

	// Evaluate curve at params and check against control points.

	for(int count = 0; count < 10; count ++)
	{
		snlPoint pt = curve.evalHmg(params[count]);

		snlVector vect = pt - points[count];

		//cout << " Point " << count << " error: " << vect.length() << "\n";

		if(maxError < vect.length()) maxError = vect.length();
	}

	cout << "Max Error = " << maxError << " - ";

	if(maxError > 1.0e-14)
		cout << "Failed\n";
	else
		cout << "Passed\n";

	delete[] points;
	delete[] params;

	return true;

}

bool testKnotRemoval()
{
	snlCurve* curve = generateCurve();

	double maxParam = curve -> maxParam();
	double minParam = curve -> minParam();
	double paramStep =(maxParam - minParam) / 10.0;

	double maxError = 0;

	for(double param = minParam + paramStep; param < maxParam; param += paramStep)
	{
		// Insert and remove knot at param up to degree times.

		int multi = curve -> knotVector().findMultiplicity(param);

		curve -> insertKnots(param, curve -> degree() - multi, true);

		for(int removalNum = 0; removalNum < curve -> degree() - multi; removalNum ++)
		{
			unsigned rSpan = curve -> knotVector().findSpan(param);

			double error = curve -> removeKnot(rSpan, 0.0);

			if(error > maxError) maxError = error;
		}

	}

	cout << "Max Error = " << maxError << " - ";

	if(maxError > 1.0e-12)
	{
		cout << "Failed\n";
		return false;
	}
	else
		cout << "Passed\n";

	delete curve;

	return true;
}

bool testDegreeElevation()
{
	cout << "\n\n";

	snlCurve* curve = generateCurve();

	bool passed = true;

	// Generate points on original curve.

	snlPoint* testPoints = new snlPoint[100];

	double maxParam = curve -> maxParam();
	double minParam = curve -> minParam();

	double paramStep =(maxParam - minParam) / 99.0;

	double param = minParam;

	for(int index = 0; index < 100; index ++)
	{
		testPoints[index] = curve -> eval(param);

		param += paramStep;
	}

	double maxError = 0.0;

	// Test elevation by n degrees.

	for(int degElev = 1; degElev <= 5; degElev ++)
	{

		snlCurve* compareCurve = new snlCurve(*curve);

		compareCurve -> elevateDegree(degElev);

		param = minParam;

		for(int index = 0; index < 100; index ++)
		{
			snlPoint testPoint = compareCurve -> eval(param);

			param += paramStep;

			double error = (testPoints[index] - testPoint).length();

			if(error > maxError) maxError = error;
		}

		cout << "\tCurve Elevation By " << degElev << " - Max Error = " << maxError << "\n";

		delete compareCurve;

		if(maxError > 1.5e-14) passed = false;
	}

	delete curve;
	delete[] testPoints;

	if(! passed)
	{
		cout << "\tFailed\n";
		return false;
	}
	else
		cout << "\tPassed\n";

	return passed;
}

bool testCurveAppend()
{
	bool passed = true;

	// Generate two curves that can be joined.

	snlCurve* curve1 = generateCurve();

	snlCurve* curve2 = generateReverseCurve();

	// Sample points along both curves to test later.

	snlPoint* testPoints1 = new snlPoint[10];

	double maxParam1 = curve1 -> maxParam();
	double minParam1 = curve1 -> minParam();

	double paramStep =(maxParam1 - minParam1) / 9.0;

	double param = minParam1;

	for(int index = 0; index < 10; index ++)
	{
		testPoints1[index] = curve1 -> eval(param);

		param += paramStep;
	}

	snlPoint* testPoints2 = new snlPoint[10];

	double maxParam2 = curve2 -> maxParam();
	double minParam2 = curve2 -> minParam();

	paramStep =(maxParam2 - minParam2) / 9.0;

	param = minParam2;

	for(int index = 0; index < 10; index ++)
	{
		testPoints2[index] = curve2 -> eval(param);

		param += paramStep;
	}

	// Append curve2 to curve1.

	int origSize = curve1 -> size();

	curve1 -> appendCurve(curve2, false);

	// Step through

	double joinParam = curve1 -> param(origSize);

	double paramStep1 = joinParam / 9.0;
	double paramStep2 =(curve1 -> maxParam() - joinParam) / 9.0;

	double param1 = curve1 -> minParam();
	double param2 = joinParam;

	double maxError = 0.0;

	for(int index = 0; index < 10; index ++)
	{
			snlPoint testPoint1 = curve1 -> eval(param1);
			snlPoint testPoint2 = curve1 -> eval(param2);

			param1 += paramStep1;
			param2 += paramStep2;

			double error1 = (testPoints1[index] - testPoint1).length();
			double error2 = (testPoints2[index] - testPoint2).length();

			if(error1 > maxError) maxError = error1;
			if(error2 > maxError) maxError = error2;
	}

	cout << "Max Error = " << maxError << " - ";

	if(maxError > 2.0e-14)
	{
		cout << "Failed\n";
		passed = false;
	}
	else
		cout << "Passed\n";

	delete curve1;
	delete curve2;

	delete[] testPoints1;
	delete[] testPoints2;

	return passed;
}

int main(int argc, char* argv[])
{
	bool allTestsPassed = true;

	cout << "\nTesting surface point interpolation ... " << flush;
	if(! test_surfaceInterp()) allTestsPassed = false;

	cout << "\nTesting surface refinement ... " << flush;
	if(! test_surfaceRefine()) allTestsPassed = false;

	cout << "\nTesting ambiguous edge detection ... " << flush;
	if(! test_ambig()) allTestsPassed = false;

	//cout << "\nTesting Surface Of Rotation ... " << flush;
	//if(! testSurfaceOfRevolution()) allTestsPassed = false;

	cout << "\nTesting Surface Point Inversion ... " << flush;
	if(! testSurfacePointInversion()) allTestsPassed = false;

	cout << "\nTesting Surface Convex Bezier Decomposition ... " << flush;
	if(! testSurfaceConvexBezierSegmentation()) allTestsPassed = false;

	cout << "\nTesting Surface Bezier Decomposition ... " << flush;
	if(! testSurfaceBezierSegmentation()) allTestsPassed = false;

	cout << "\nTesting Surface Knot Insertion ... " << flush;
	if(! testSurfaceKnotInsert()) allTestsPassed = false;

	cout << "\nTesting Derivative Evaluation ... " << flush;
	if(! testDerivEval()) allTestsPassed = false;

	cout << "\nTesting Surface Projection ... " << flush;
	if(! testSurfaceProjection()) allTestsPassed = false;

	cout << "\nTesting snlSquareLinear ... " << flush;
	if(testSquareLinear())
		cout << "Passed\n";
	else
	{
		cout << "Failed\n";
		allTestsPassed = false;
	}

	cout << "\nTesting Curve Interpolation\n";
	if(! testCurveInterpolation()) allTestsPassed = false;

	cout << "\nTesting Curve Knot Removal ... " << flush;
	if(! testKnotRemoval()) allTestsPassed = false;

	cout << "\nTesting Degree Elevation ... " << flush;
	if(! testDegreeElevation()) allTestsPassed = false;

	cout << "\nTesting Curve Append ... " << flush;
	if(! testCurveAppend()) allTestsPassed = false;

	cout << "\n";

	if(allTestsPassed)
		cout << "All tests have passed :-)\n\n";
	else
		cout << "*** Some tests have not passed ***\n\n";
}


