EXEC+=mc
EXEC+=check
EXEC+=min
EXEC+=study
EXEC+=load

mc_SRCS=    mc.cpp    MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareACSL.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp
min_SRCS=   min.cpp   MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareACSL.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp MCSim.cpp VMCMinimization.cpp Interpolation.cpp VMCInterpolation.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp 
check_SRCS= check.cpp MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareACSL.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp
study_SRCS= study.cpp MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareACSL.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp Directory.cpp Analyse.cpp AnalyseEnergy.cpp AnalyseChain.cpp AnalyseMagnetization.cpp AnalyseHoneycomb.cpp AnalyseLadder.cpp MCSim.cpp VMCMinimization.cpp Interpolation.cpp
load_SRCS=  load.cpp  MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareACSL.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp Directory.cpp Analyse.cpp AnalyseEnergy.cpp AnalyseChain.cpp AnalyseMagnetization.cpp AnalyseHoneycomb.cpp AnalyseLadder.cpp MCSim.cpp VMCMinimization.cpp Interpolation.cpp

#-----------------------------------------------------------------

CXX = g++ -std=c++11

ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -fopenmp
FIT = /home/jdufour/travail/cpp-dev/fit/

BUILD=build

CXXFLAGS = -I$(FIT) $(ERRORS) $(OPTION)
LDFLAGS  = $(FIT)libcminpack.a $(LAPACK) $(ERRORS) $(OPTION)

SRCS=$(wildcard *.cpp)

all: OPTION += -O3 -DNDEBUG
all:$(EXEC)
	cp mc ../sim
	cp check ../sim
	cp min ../sim
	cp load ../sim
	cp study ..

debug: OPTION += -ggdb 
debug:$(EXEC)
	cp mc ../sim
	cp check ../sim
	cp min ../sim
	cp load ../sim
	cp study ..

.SECONDEXPANSION:
$(EXEC): $$(patsubst %.cpp, $(BUILD)/%.o, $$($$@_SRCS)) 
	@echo Links $(notdir $^)
	$(CXX) -o $@ $^ $(LDFLAGS) 

$(BUILD)/%.o:%.cpp
	@echo Creates $(notdir $@)
	$(CXX) -MD -c $(CXXFLAGS) $< -o $@

-include $(addprefix $(BUILD)/,$(SRCS:.cpp=.d))

clean:
	rm -f $(BUILD)/* $(EXEC)

ref:
	@echo Create the documentation
	@doxygen doxygen/Doxyfile > doxygen/log
	@firefox doxygen/html/files.html &
