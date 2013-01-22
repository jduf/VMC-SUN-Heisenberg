CXX = g++ 
CC = $(CXX) 

DEBUG = -g 
ERRORS = -pedantic -Wall -Wextra
LAPACK = -llapack -lblas
OPT = -O3
#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2 
CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG)
#CXXFLAGS = -O3 -fopenmp -msse2

#LIBFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS= $(LAPACK) $(DEBUG)
#LIBFLAGS = -O3 -fopenmp -msse2

all:test

test:test.o State.o System.o Chrono.o Vecteur.o Matrice.o Lapack.o
	$(CXX) -o $@ $^ $(LDFLAGS)

State.o:State.cpp State.hpp System.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^

System.o:System.cpp System.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^

Chrono.o:Chrono.cpp Chrono.hpp
	$(CXX) -c $(CXXFLAGS) $^

Vecteur.o:Vecteur.cpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

Matrice.o:Matrice.cpp Matrice.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^
	
Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

test.o:test.cpp Save.hpp State.hpp System.hpp Chrono.hpp Matrice.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

clean:
	rm *.o *.gch test
