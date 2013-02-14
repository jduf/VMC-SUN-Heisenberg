CXX = g++ 
CC = $(CXX) 

#DEBUG = -DNDEBUG
DEBUG = -pg
ERRORS = -Wall -Wextra
LAPACK = -llapack -lblas -pedantic
OPTION = 

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS= $(LAPACK) $(DEBUG) $(OPTION)

all:test

test:test.o System.o Lapack.o Read.o
	$(CXX) -o $@ $^ $(LDFLAGS)

Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp 
	$(CXX) -c $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Matrice.hpp 
	$(CXX) -c $(CXXFLAGS) $^

System.o:System.cpp System.hpp Vecteur.hpp Matrice.hpp Lapack.hpp Read.hpp
	$(CXX) -c $(CXXFLAGS) $^

test.o:test.cpp Save.hpp Chrono.hpp Parseur.hpp System.hpp 
	$(CXX) -c $(CXXFLAGS) $^ 

clean:
	rm *.o *.gch test
