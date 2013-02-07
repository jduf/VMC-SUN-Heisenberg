CXX = g++ 
CC = $(CXX) 

DEBUG =
ERRORS = -Wall -Wextra
LAPACK = -llapack -lblas
OPTION = -O3

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS= $(LAPACK) $(DEBUG) $(OPTION)

all:test

test:test.o Lapack.o Matrice.o System.o
	$(CXX) -o $@ $^ $(LDFLAGS)

Matrice.o:Matrice.cpp Matrice.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

System.o:System.cpp System.hpp Matrice.hpp Lapack.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp 
	$(CXX) -c $(CXXFLAGS) $^

test.o:test.cpp Save.hpp Chrono.hpp System.hpp 
	$(CXX) -c $(CXXFLAGS) $^ 

clean:
	rm *.o *.gch test
