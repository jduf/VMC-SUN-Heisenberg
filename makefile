CXX = g++ 
CC = $(CXX) 

DEBUG = -DNDEBUG
#DEBUG = -pg
ERRORS = -Wall -Wextra
LAPACK = -llapack -lblas -pedantic
OPTION = -O3

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS= $(LAPACK) $(DEBUG) $(OPTION)

all:mc setup
########
# mc
#######
mc:MC.o Lapack.o Read.o Write.o
	$(CXX) -o $@ $^ $(LDFLAGS)

Matrice.o:Matrice.cpp Matrice.hpp Vecteur.hpp
	$(CXX) -c $(CXXFLAGS) $^

Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp 
	$(CXX) -c $(CXXFLAGS) $^

MC.o:MC.cpp Chrono.hpp Parseur.hpp System.hpp Read.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

########
# setup
########
setup:Setup.o Write.o Read.o
	$(CXX) -o $@ $^ $(LDFLAGS)

Setup.o:Setup.cpp Read.hpp Write.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

########
# commun
########
Write.o:Write.cpp Write.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^

########
# divers
########
clean:
	rm *.o *.gch mc setup
