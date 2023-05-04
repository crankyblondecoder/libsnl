# libSNL Test Directory Makefile
# ------------------------------

VPATH = ../../src/

cname = g++
cflags = -Wall -fPIC -g

CXX = $(cname)
CXXFLAGS = $(cflags)

test : snlTest.cpp
	$(cname) $(cflags) snlTest.cpp -o snlTest $(mainObjFiles)
