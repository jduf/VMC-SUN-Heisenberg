CXX = g++

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

EXEC=mc min analyse check
SRCS=$(wildcard *.cpp)
BUILD=build

all:$(EXEC)
	cp $(EXEC) ../sim

mc:mc.o CreateSystem.o ChainFermi.o ChainPolymerized.o SquarePiFlux.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o
	@echo Links $@
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

min:min.o CreateSystem.o ChainFermi.o ChainPolymerized.o SquarePiFlux.o Minimization.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o
	@echo Links $@
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

analyse:analyse.o Parseur.o Read.o Write.o Header.o RST.o Gnuplot.o Directory.o
	@echo Links $@
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

check:check.o CreateSystem.o ChainFermi.o ChainPolymerized.o SquarePiFlux.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o
	@echo Links $@
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

#pso:pso.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o  PSTricks.o PSO.o PSOFermionic.o CreateSystem.o ChainFermi.o ChainDimerized.o TriangleJastrow.o SquareJastrow.o SquareSU2PhiFlux.o SquarePiFlux.o SquareMu.o SquareFermi.o HoneycombSU3.o HoneycombSU4.o Write.o Read.o TriangleFermi.o TriangleMu.o TrianglePhi.o  ChainPolymerized.o
#	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)
#psomc:psomc.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o SquareJastrow.o TriangleJastrow.o  PSTricks.o PSO.o PSOMonteCarlo.o
#	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)
#mcnu:mcnu.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o  SquareJastrow.o TriangleJastrow.o PSTricks.o
#	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

%.o:%.cpp
	@echo Creates $@
	$(CXX) -MD -c $(CXXFLAGS) $(NOASSERT)  $< -o $@
	@mv $(<:.cpp=.d) $(BUILD)

-include $(addprefix $(BUILD)/,$(SRCS:.cpp=.d))

clean:
	rm -f *.gch *.o $(BUILD)/* $(EXEC)

ref:
	@echo Create the documentation
	@doxygen doxygen/Doxyfile
	@firefox doxygen/html/files.html &
