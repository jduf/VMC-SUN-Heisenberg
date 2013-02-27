CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS = -pedantic -Wall -Wextra 
#LDLIBS = -O3 -fopenmp -msse2

all:test

test.o:test.cpp Read.hpp Write.hpp Matrice.hpp Vecteur.hpp Array2D.hpp
Read.o:Read.cpp Read.hpp Matrice.hpp Array2D.hpp
Write.o:Write.cpp Write.hpp Matrice.hpp Array2D.hpp
test:test.o Read.o Write.o

ref:
	doxygen Doxyfile

clean:
	rm *.o test
