# libSNL Root Makefile
# --------------------

export objFiles := snlCircularOffsetCurve.o \
       snlCtrlPoint.o snlCtrlPointNet.o snlCtrlPointNetCurve.o snlCtrlPointNetSurface.o snlCurve.o \
       snlKnotVector.o \
       snlMatrix_4x4.o snlMeshable.o \
       snlNurbsCommon.o \
       snlPoint.o \
       snlSquareLinear.o snlSurface.o snlSurface_pointLoop.o \
       snlSurface_projection.o snlSurfaceOfRevolution.o \
       snlTransform.o snlTriangleMesh.o \
       snlUtil.o \
       snlVector.o snlVertex.o snlVertexNet.o

.PHONY : debug
debug :
	$(MAKE) -C debug debug

.PHONY : clean
clean :
	- $(MAKE) -C debug clean
	- $(MAKE) -C test clean
	- $(MAKE) -C release clean

.PHONY : test
test : debug
	$(MAKE) -C test test

.PHONY : release
release :
	$(MAKE) -C release release