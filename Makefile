#/bin/sh
LIB_CGNS=#-lcgns 
CGNS_INCLUDE=#/usr/local/include

SRC_FILES = src/basis.f90 \
	src/evaluate.f90 \
	src/knots.f90 \
	src/jacobian.f90 \
	src/paramuni.f90

default:	
	f2py --fcompiler=gnu95 -c ${SRC_FILES} -m MBIlib
	mv MBIlib.so MBI/
