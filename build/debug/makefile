# libSNL Source Directory Debug Makefile
# -------------------------------------

# Tell make to search other paths. Separate multiple with ":"
VPATH = ../../src/

libName = ../libSNL_debug.so.0.2

cname = g++
cflags = -Wall -fPIC -g

CXX = $(cname)
CXXFLAGS = $(cflags)

.PHONY : debug
debug : $(objFiles)
	@echo
	@echo "*** Building - Debug ***"
	@echo
	$(cname) $(cflags) -shared -o $(libName) $(objFiles)

# GNU Make will automatically update included files.
ifneq ($(MAKECMDGOALS),clean)
    include $(objFiles:.o=.dep)
endif

%.dep : %.cpp %.h
	$(cname) -MM $< > $@

.PHONY : clean
clean :
	rm *.o *.dep $(libName)
