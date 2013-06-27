CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDLIBS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test2

Lapack.o:Lapack.cpp Lapack.hpp Matrix.hpp Vecteur.hpp
test.o:test.cpp Vecteur.hpp Matrice.hpp Lapack.hpp 
test2.o:test2.cpp Matrix.hpp Vecteur.hpp

test:test.o Lapack.o

test2:test2.o Lapack.o

ref:
	doxygen Doxyfile

clean:
	rm *.o test
