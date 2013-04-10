CXX = g++ 
CC = $(CXX) 

#DEBUG = -pg
DEBUG = -DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS  = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

all:mc cs check
	cp mc cs check sim/

#############
# monte-carlo
#############
mc:mc.o Lapack.o Rand.o Read.o Write.o Header.o RST.o Chrono.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

Chrono.o:Chrono.cpp Chrono.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

mc.o:mc.cpp MonteCarlo.hpp System.hpp Lapack.hpp Read.hpp Array2D.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

########
# create
########
cs:cs.o CreateSystem.o Lapack.o Read.o Write.o Header.o RST.o 
	$(CXX) -o $@ $^ $(LAPACK) $(OPTION)

CreateSystem.o:CreateSystem.cpp CreateSystem.hpp Read.hpp Write.hpp Array2D.hpp Matrice.hpp Lapack.hpp RST.hpp 
	$(CXX) -c $(LAPACK) $(OPTION) $^ 

cs.o:cs.cpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(LAPACK) $(OPTION) $^ 

#######
# check
#######

check:check.o Read.o Header.o RST.o Write.o
	$(CXX) -o $@ $^ $(OPTION)

check.o:check.cpp Matrice.hpp Array2D.hpp Read.hpp
	$(CXX) -c $(OPTION) $^ 
	
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
	rm *.o *.gch mc cs check

ref:
	doxygen Doxyfile
