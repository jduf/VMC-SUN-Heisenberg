CXX = g++ 
CC = $(CXX) 

DEBUG =  
ERRORS = -Wall -Wextra
LAPACK = -llapack -lblas
OPTION = 

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS= $(LAPACK) $(DEBUG) $(OPTION)

all:test

test:test.o State.o System.o Vecteur.o Matrice.o Lapack.o
	$(CXX) -o $@ $^ $(LDFLAGS)

State.o:State.cpp State.hpp System.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^

System.o:System.cpp System.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^

Vecteur.o:Vecteur.cpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

Matrice.o:Matrice.cpp Matrice.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^
	
Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

test.o:test.cpp Save.hpp Chrono.hpp State.hpp System.hpp Chrono.hpp Matrice.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

clean:
	rm *.o *.gch test
