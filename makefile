CXX = g++ 
CC = $(CXX) 

DEBUG = -DNDEBUG
#DEBUG = -pg
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas -pedantic
OPTION = -O3

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS  = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

all:mc cs
########
# mc
#######
mc:mc.o MonteCarlo.o Lapack.o Rand.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

MonteCarlo.o:MonteCarlo.cpp MonteCarlo.hpp Chrono.hpp System.hpp Read.hpp Array2D.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

mc.o:mc.cpp MonteCarlo.hpp Chrono.hpp System.hpp Read.hpp Array2D.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

########
# create
########
cs:cs.o CreateSystem.o Lapack.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LAPACK) $(OPTION)

CreateSystem.o:CreateSystem.cpp CreateSystem.hpp Read.hpp Write.hpp Array2D.hpp Matrice.hpp Lapack.hpp
	$(CXX) -c $(LAPACK) $(OPTION) $^ 

cs.o:cs.cpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(LAPACK) $(OPTION) $^ 

########
# commun
########
Write.o:Write.cpp Write.hpp Header.hpp Matrice.hpp Array2D.hpp
	$(CXX) -c $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Header.hpp Matrice.hpp Array2D.hpp
	$(CXX) -c $(CXXFLAGS) $^
	
Header.o:Header.cpp Header.hpp Matrice.hpp Array2D.hpp
	$(CXX) -c $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp
	$(CXX) -c $(CXXFLAGS) $^

Lapack.o:Lapack.cpp Lapack.hpp Matrice.hpp Vecteur.hpp 
	$(CXX) -c $(CXXFLAGS) $^
########
# divers
########
clean:
	rm *.o *.gch mc cs

ref:
	doxygen Doxyfile
