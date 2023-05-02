# libSNL Test Directory Makefile
# ------------------------------

mainObjFiles = $(addprefix ../, $(objFiles))

CXX = $(cname)
CXXFLAGS = $(cflags)

test : snlTest.cpp
	$(cname) $(cflags) snlTest.cpp -o snlTest $(mainObjFiles)
