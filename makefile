CXX = g++ -g
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
mc:mc.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o 
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

mc.o:mc.cpp Parseur.hpp MonteCarlo.hpp System.hpp Read.hpp Matrix.hpp Lapack.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^ 

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^ 

########
# create
########
cs:cs.o Chain.o SquareSU2PhiFlux.o SquarePiFlux.o SquareMu.o SquareFermi.o HoneycombSU4.o Parseur.o Write.o Read.o Lapack.o RST.o Header.o Gnuplot.o PSTricks.o TriangleFermi.o TriangleMu.o TrianglePhi.o
	$(CXX) -o $@ $^ $(LDFLAGS)

cs.o:cs.cpp Parseur.hpp Chain.hpp Square.hpp Honeycomb.hpp Triangle.hpp
	$(CXX) -c $(CXXFLAGS) $^

Chain.o:Chain.cpp Chain.hpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquarePiFlux.o:SquarePiFlux.cpp SquarePiFlux.hpp Square.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^
	
SquareSU2Phiflux.o:SquareSU2PhiFlux.cpp SquareSU2PhiFlux.hpp Square.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareMu.o:SquareMu.cpp SquareMu.hpp Square.hpp CreateSystem.hpp Parseur.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareFermi.o:SquareFermi.cpp SquareFermi.hpp Square.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^

TriangleFermi.o:TriangleFermi.cpp TriangleFermi.hpp Triangle.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^

TriangleMu.o:TriangleMu.cpp TriangleMu.hpp Triangle.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^

TrianglePhi.o:TrianglePhi.cpp TrianglePhi.hpp Triangle.hpp CreateSystem.hpp Parseur.hpp Gnuplot.hpp Lapack.hpp PSTricks.hpp
	$(CXX) -c $(CXXFLAGS) $^

HoneycombSU4.o:HoneycombSU4.cpp HoneycombSU4.hpp Honeycomb.hpp CreateSystem.hpp Parseur.hpp
	$(CXX) -c $(CXXFLAGS) $^

Gnuplot.o:Gnuplot.cpp Gnuplot.hpp Write.hpp RST.hpp Header.hpp Time.hpp Matrix.hpp Vector.hpp
	$(CXX) -c $(CXXFLAGS) $^
	
PSTricks.o:PSTricks.cpp PSTricks.hpp Write.hpp Linux.hpp Vector.hpp
	$(CXX) -c $(CXXFLAGS) $^

#######
# check
#######
check:check.o Read.o Header.o RST.o Write.o
	$(CXX) -o $@ $^ $(LDFLAGS)

check.o:check.cpp Read.hpp Matrix.hpp Container.hpp
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

Lapack.o:Lapack.cpp Lapack.hpp Matrix.hpp Vector.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

########
# divers
########
clean:
	rm *.o *.gch mc cs check

ref:
	doxygen doxygen/Doxyfile
	firefox doxygen/html/files.html
