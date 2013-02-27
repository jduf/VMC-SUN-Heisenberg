CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDLIBS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test 

Lapack.o:Lapack.cpp Lapack.hpp
test.o:test.cpp Vecteur.hpp Matrice.hpp Lapack.hpp 

test:test.o Lapack.o

ref:
	doxygen Doxyfile

clean:
	rm *.o test
