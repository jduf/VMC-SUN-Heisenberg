#EXEC+=mc
EXEC+=check
#EXEC+=study
#EXEC+=psomc

mc_SRCS=    mc.cpp    System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp TriangleFermi.cpp SquareFermi.cpp SquarePiFlux.cpp SquareJastrow.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp
check_SRCS= check.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp TriangleFermi.cpp SquareFermi.cpp SquarePiFlux.cpp SquareJastrow.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp
study_SRCS= study.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp TriangleFermi.cpp SquareFermi.cpp SquarePiFlux.cpp SquareJastrow.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp Directory.cpp Analyse.cpp AnalyseEnergy.cpp AnalyseChain.cpp AnalyseMagnetization.cpp AnalyseHoneycomb.cpp
psomc_SRCS= psomc.cpp System.cpp IOSystem.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp TriangleFermi.cpp SquareFermi.cpp SquarePiFlux.cpp SquareJastrow.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Honeycomb0pp.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Fit.cpp Gnuplot.cpp PSO.cpp PSOFermionic.cpp

#-----------------------------------------------------------------

CXX = g++ -std=c++11

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 -fopenmp

BUILD=build

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

SRCS=$(wildcard *.cpp)

all:$(EXEC)
	cp mc ../sim
	cp check ../sim
	cp psomc ../sim
	cp study ..

.SECONDEXPANSION:
$(EXEC): $$(patsubst %.cpp, $(BUILD)/%.o, $$($$@_SRCS)) 
	@echo Links $(notdir $^)
	$(CXX) -o $@ $^ libcminpack.a $(LDFLAGS) $(NOASSERT)

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
