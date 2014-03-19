EXEC=mc check study

mc_SRCS=   mc.cpp    CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp Lapack.cpp Parseur.cpp Read.cpp Write.cpp Header.cpp RST.cpp PSTricks.cpp Rand.cpp
#min_SRCS= min.cpp   CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp Lapack.cpp Parseur.cpp Read.cpp Write.cpp Header.cpp RST.cpp PSTricks.cpp Rand.cpp  Minimization.cpp
check_SRCS=check.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp Lapack.cpp Parseur.cpp Read.cpp Write.cpp Header.cpp RST.cpp PSTricks.cpp
study_SRCS=study.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp Lapack.cpp Parseur.cpp Read.cpp Write.cpp Header.cpp RST.cpp PSTricks.cpp RSTFile.cpp Gnuplot.cpp Directory.cpp

#pso:pso.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o  PSTricks.o PSO.o PSOFermionic.o CreateSystem.o ChainFermi.o ChainDimerized.o TriangleJastrow.o SquareJastrow.o SquareSU2PhiFlux.o SquarePiFlux.o SquareMu.o SquareFermi.o HoneycombSU3.o HoneycombSU4.o Write.o Read.o TriangleFermi.o TriangleMu.o TrianglePhi.o  ChainPolymerized.o
#	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)
#psomc:psomc.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o SquareJastrow.o TriangleJastrow.o  PSTricks.o PSO.o PSOMonteCarlo.o
#	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)
#mcnu:mcnu.o Parseur.o Lapack.o Rand.o Read.o Write.o Header.o RST.o  SquareJastrow.o TriangleJastrow.o PSTricks.o
#	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

#-----------------------------------------------------------------

CXX = g++

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

BUILD=build

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

SRCS=$(wildcard *.cpp)

all:$(EXEC)
	cp mc check ../sim
	cp study ..


.SECONDEXPANSION:
$(EXEC): $$(patsubst %.cpp, %.o, $$($$@_SRCS)) 
	@echo Links $^
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

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
