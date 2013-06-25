CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDLIBS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test test2

Lapack.o:Lapack.cpp Lapack.hpp
test.o:test.cpp Vecteur.hpp Matrice.hpp Lapack.hpp 
test2.o:test2.cpp Array.hpp Vecteur.hpp 

test:test.o Lapack.o

test2:test2.o 

ref:
	doxygen Doxyfile

clean:
	rm *.o test
