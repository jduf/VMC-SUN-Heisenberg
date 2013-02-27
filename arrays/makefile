CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDLIBS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test 


test.o:test.cpp Array1D.hpp Array2D.hpp Read.hpp Write.hpp
Read.o:Read.cpp Read.hpp Matrice.hpp
Write.o:Write.cpp Write.hpp Matrice.hpp

test:test.o Read.o Write.o

ref:
	doxygen Doxyfile

clean:
	rm *.o test
