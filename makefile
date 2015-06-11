EXEC+=mc
EXEC+=check
EXEC+=study
EXEC+=psomc
EXEC+=psplinemc
EXEC+=min

mc_SRCS=    mc.cpp            MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerizedJJp.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp
check_SRCS= check.cpp         MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerizedJJp.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp
study_SRCS= study.cpp         MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerizedJJp.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp Directory.cpp Analyse.cpp AnalyseEnergy.cpp AnalyseChain.cpp AnalyseMagnetization.cpp AnalyseHoneycomb.cpp
psomc_SRCS= psomc.cpp         MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerizedJJp.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp MCSim.cpp VMCMinimization.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp 
min_SRCS=   min.cpp           MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerizedJJp.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp MCSim.cpp VMCMinimization.cpp PSpline.cpp VMCSpline.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp 
psplinemc_SRCS= psplinemc.cpp MonteCarlo.cpp MCSystem.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerizedJJp.cpp ChainPolymerized.cpp LadderFermi.cpp LadderFree.cpp SquareFermi.cpp SquarePiFlux.cpp SquareFreeComplex.cpp SquareJastrow.cpp TriangleFermi.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp MCSim.cpp VMCMinimization.cpp PSpline.cpp VMCSpline.cpp

#-----------------------------------------------------------------

CXX = g++ -std=c++11

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp
FIT = /home/jdufour/travail/cpp-dev/fit/

BUILD=build

CXXFLAGS = -I$(FIT) $(ERRORS) $(OPTION)
LDFLAGS  = $(FIT)libcminpack.a $(LAPACK) $(ERRORS) $(OPTION)

SRCS=$(wildcard *.cpp)

all:$(EXEC)
	cp mc ../sim
	cp check ../sim
	cp psomc ../sim
	cp min ../sim
	cp study ..

.SECONDEXPANSION:
$(EXEC): $$(patsubst %.cpp, $(BUILD)/%.o, $$($$@_SRCS)) 
	@echo Links $(notdir $^)
	$(CXX) -o $@ $^ $(LDFLAGS) $(NOASSERT)

$(BUILD)/%.o:%.cpp
	@echo Creates $(notdir $@)
	$(CXX) -MD -c $(CXXFLAGS) $(NOASSERT)  $< -o $@

-include $(addprefix $(BUILD)/,$(SRCS:.cpp=.d))

clean:
	rm -f $(BUILD)/* $(EXEC)

ref:
	@echo Create the documentation
	@doxygen doxygen/Doxyfile
	@firefox doxygen/html/files.html &
