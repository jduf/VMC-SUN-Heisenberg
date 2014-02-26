CXX = g++ -g
CC = $(CXX) 

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

CXXFLAGS = $(LAPACK) $(ERRORS) $(OPTION)

LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

all:mc min
	cp mc min ../sim

#############
# monte-carlo
#############
mc:mc.o CreateSystem.o ChainFermi.o ChainPolymerized.o SquarePiFlux.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

mc.o:mc.cpp CreateSystem.hpp MonteCarlo.hpp System.hpp SystemBosonic.hpp SystemFermionic.hpp Chain.hpp ChainPolymerized.hpp ChainFermi.hpp Square.hpp SquarePiFlux.hpp Parseur.hpp Vector.hpp Matrix.hpp Lapack.hpp Read.hpp Write.hpp RST.hpp Header.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

##############
# minimization
##############
min:min.o CreateSystem.o ChainFermi.o ChainPolymerized.o SquarePiFlux.o Minimization.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

min.o:min.cpp Minimization.hpp CreateSystem.hpp MonteCarlo.hpp System.hpp SystemBosonic.hpp SystemFermionic.hpp Chain.hpp ChainPolymerized.hpp ChainFermi.hpp Square.hpp SquarePiFlux.hpp Parseur.hpp Vector.hpp Matrix.hpp Lapack.hpp Read.hpp Write.hpp Write.hpp RST.hpp Header.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

Minimization.o:Minimization.cpp Minimization.hpp Parseur.hpp MonteCarlo.hpp System.hpp SystemFermionic.hpp SystemBosonic.hpp Read.hpp  Matrix.hpp Lapack.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

##############
# systems
##############
CreateSystem.o:CreateSystem.cpp CreateSystem.hpp GenericSystem.hpp Parseur.hpp Chain.hpp ChainFermi.hpp ChainPolymerized.hpp 
	$(CXX) -c $(CXXFLAGS) $^

ChainFermi.o:ChainFermi.cpp ChainFermi.hpp Chain.hpp GenericSystem.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

ChainDimerized.o:ChainDimerized.cpp ChainDimerized.hpp Chain.hpp GenericSystem.hpp  Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

ChainPolymerized.o:ChainPolymerized.cpp ChainPolymerized.hpp Chain.hpp GenericSystem.hpp  Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquarePiFlux.o:SquarePiFlux.cpp SquarePiFlux.hpp Square.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareSU2Phiflux.o:SquareSU2PhiFlux.cpp SquareSU2PhiFlux.hpp Square.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareJastrow.o:SquareJastrow.cpp SquareJastrow.hpp Square.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

TriangleJastrow.o:TriangleJastrow.cpp TriangleJastrow.hpp Triangle.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareMu.o:SquareMu.cpp SquareMu.hpp Square.hpp GenericSystem.hpp   Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

SquareFermi.o:SquareFermi.cpp SquareFermi.hpp Square.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

TriangleFermi.o:TriangleFermi.cpp TriangleFermi.hpp Triangle.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

TriangleMu.o:TriangleMu.cpp TriangleMu.hpp Triangle.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

TrianglePhi.o:TrianglePhi.cpp TrianglePhi.hpp Triangle.hpp GenericSystem.hpp  Gnuplot.hpp Lapack.hpp PSTricks.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

HoneycombSU4.o:HoneycombSU4.cpp HoneycombSU4.hpp Honeycomb.hpp GenericSystem.hpp  Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

HoneycombSU3.o:HoneycombSU3.cpp HoneycombSU3.hpp Honeycomb.hpp GenericSystem.hpp  Container.hpp
	$(CXX) -c $(CXXFLAGS) $^

########
# commun
########
Gnuplot.o:Gnuplot.cpp Gnuplot.hpp Write.hpp RST.hpp Header.hpp Time.hpp Matrix.hpp Vector.hpp
	$(CXX) -c $(CXXFLAGS) $^

PSTricks.o:PSTricks.cpp PSTricks.hpp Write.hpp Linux.hpp Vector.hpp
	$(CXX) -c $(CXXFLAGS) $^

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

#############################
# particle-swarm optimisation For fermionic thing
#############################
pso:pso.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o Container.o PSTricks.o PSO.o PSOFermionic.o CreateSystem.o ChainFermi.o ChainDimerized.o TriangleJastrow.o SquareJastrow.o SquareSU2PhiFlux.o SquarePiFlux.o SquareMu.o SquareFermi.o HoneycombSU3.o HoneycombSU4.o Write.o Read.o TriangleFermi.o TriangleMu.o TrianglePhi.o Container.o ChainPolymerized.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

pso.o:pso.cpp Parseur.hpp MonteCarlo.hpp System.hpp SystemFermionic.hpp SystemBosonic.hpp Read.hpp  Matrix.hpp Lapack.hpp Container.hpp PSO.hpp PSOFermionic.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

PSOFermionic.o:PSOFermionic.cpp PSOFermionic.hpp PSO.hpp MonteCarlo.hpp CreateSystem.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

#############################
# particle-swarm optimisation
#############################
psomc:psomc.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o SquareJastrow.o TriangleJastrow.o Container.o PSTricks.o PSO.o PSOMonteCarlo.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

psomc.o:psomc.cpp Parseur.hpp MonteCarlo.hpp System.hpp SystemFermionic.hpp SystemBosonic.hpp Read.hpp  Matrix.hpp Lapack.hpp Container.hpp PSO.hpp PSOMonteCarlo.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

PSO.o:PSO.cpp PSO.hpp Rand.hpp Write.hpp Read.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

PSOMonteCarlo.o:PSOMonteCarlo.cpp PSOMonteCarlo.hpp PSO.hpp MonteCarlo.hpp SquareJastrow.hpp TriangleJastrow.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

#######
# mc nu
#######
mcnu:mcnu.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o Container.o SquareJastrow.o TriangleJastrow.o PSTricks.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

mcnu.o:mcnu.cpp Parseur.hpp MonteCarlo.hpp System.hpp SystemFermionic.hpp SystemBosonic.hpp Read.hpp  Matrix.hpp Lapack.hpp Container.hpp PSO.hpp PSOMonteCarlo.hpp TriangleJastrow.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

########
# useless
########
ExtractSystem.o:ExtractSystem.cpp ExtractSystem.hpp Container.hpp
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

Container.o:Container.cpp Container.hpp 
	$(CXX) -c $(CXXFLAGS) $(NOASSERT) $^

check:check.o ExtractSystem.o Read.o Header.o RST.o Write.o Container.o
	$(CXX) -o $@ $^ $(LDFLAGS)

check.o:check.cpp ExtractSystem.hpp
	$(CXX) -c $(CXXFLAGS) $^

clean:
	rm *.o *.gch mc cs check psomc mcnu

ref:
	doxygen doxygen/Doxyfile
	firefox doxygen/html/files.html &
