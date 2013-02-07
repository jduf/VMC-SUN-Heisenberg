CXX = g++ 
CC = $(CXX) 

DEBUG =
ERRORS = -Wall -Wextra
LAPACK = -llapack -lblas
OPTION = -O3

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS= $(LAPACK) $(DEBUG) $(OPTION)

all:test

test:test.o System.o Lapack.o
	$(CXX) -o $@ $^ $(LDFLAGS)

Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp 
	$(CXX) -c $(CXXFLAGS) $^

System.o:System.cpp System.hpp Vecteur.hpp Matrice.hpp Lapack.hpp
	$(CXX) -c $(CXXFLAGS) $^

test.o:test.cpp Save.hpp Chrono.hpp Parseur.hpp System.hpp 
	$(CXX) -c $(CXXFLAGS) $^ 

clean:
	rm *.o *.gch test
