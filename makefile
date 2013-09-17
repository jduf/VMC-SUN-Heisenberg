CXX = g++ 
CC = $(CXX) 

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

CXXFLAGS = $(LAPACK) $(ERRORS) $(OPTION)

LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

all:mc cs check
	cp mc cs check sim

#############
# monte-carlo
#############
mc:mc.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o Chrono.o 
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

mc.o:mc.cpp Parseur.hpp MonteCarlo.hpp System.hpp Read.hpp Matrix.hpp Lapack.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^ 

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^ 

Chrono.o:Chrono.cpp Chrono.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^ 
	
########
# create
########
cs:cs.o Chain.o SquarePiFlux.o SquareMu.o SquareFermi.o HoneycombSU4.o Parseur.o Write.o Read.o Lapack.o RST.o Header.o Gnuplot.o
	$(CXX) -o $@ $^ $(LDFLAGS)

cs.o:cs.cpp Parseur.hpp Chain.hpp Square.hpp Honeycomb.hpp 
	$(CXX) -c $(CXXFLAGS) $^

Chain.o:Chain.cpp Chain.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquarePiFlux.o:SquarePiFlux.cpp SquarePiFlux.hpp Square.hpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(CXXFLAGS) $^
	
SquareMu.o:SquareMu.cpp SquareMu.hpp Square.hpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareFermi.o:SquareFermi.cpp SquareFermi.hpp Square.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp
	$(CXX) -c $(CXXFLAGS) $^

HoneycombSU4.o:HoneycombSU4.cpp HoneycombSU4.hpp Honeycomb.hpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(CXXFLAGS) $^

Gnuplot.o:Gnuplot.cpp Gnuplot.hpp Write.hpp RST.hpp Header.hpp Time.hpp Matrix.hpp 
	$(CXX) -c  $(CXXFLAGS) $^
	
#######
# check
#######
check:check.o Read.o Header.o RST.o Write.o
	$(CXX) -o $@ $^ $(LDFLAGS)

check.o:check.cpp Read.hpp Matrix.hpp 
	$(CXX) -c $(CXXFLAGS) $^
	
########
# commun
########
Parseur.o:Parseur.cpp Parseur.hpp 
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

Write.o:Write.cpp Write.hpp Header.hpp Matrix.hpp 
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

Read.o:Read.cpp Read.hpp Header.hpp Matrix.hpp 
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^
	
Header.o:Header.cpp Header.hpp Matrix.hpp Time.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

RST.o:RST.cpp RST.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

Lapack.o:Lapack.cpp Lapack.hpp Matrix.hpp  
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

########
# divers
########
clean:
	rm *.o *.gch mc cs check

ref:
	doxygen Doxyfile
	firefox html/files.html
