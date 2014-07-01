EXEC=mc check study

mc_SRCS=   mc.cpp    System.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp PSTricks.cpp Rand.cpp Gnuplot.cpp
check_SRCS=check.cpp System.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp PSTricks.cpp Rand.cpp Gnuplot.cpp
study_SRCS=study.cpp System.cpp CreateSystem.cpp ChainFermi.cpp ChainPolymerized.cpp SquarePiFlux.cpp KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp Lapack.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp PSTricks.cpp Analyse.cpp AnalyseParameter.cpp AnalyseMagnetization.cpp RSTFile.cpp Gnuplot.cpp Directory.cpp

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
	cp study ..
	cp mc check ../sim


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
