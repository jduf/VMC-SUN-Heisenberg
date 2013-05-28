CXX = g++ 
CC = $(CXX) 

#DEBUG = -pg
DEBUG = -DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

CXXFLAGS = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

LDFLAGS  = $(LAPACK) $(ERRORS) $(DEBUG) $(OPTION)

all:mc
	cp mc ../sim

#############
# monte-carlo
#############
mc:mc.o Lapack.o Rand.o Read.o Write.o Header.o RST.o Chrono.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

mc.o:mc.cpp MonteCarlo.hpp System.hpp Lapack.hpp Read.hpp Array2D.hpp Matrice.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

Chrono.o:Chrono.cpp Chrono.hpp
	$(CXX) -c $(CXXFLAGS) $^ 
	
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
	rm *.o *.gch mc 

ref:
	doxygen Doxyfile
