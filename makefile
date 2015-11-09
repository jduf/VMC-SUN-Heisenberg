EXEC+=mc
EXEC+=check
EXEC+=min
EXEC+=study
EXEC+=load

CHAIN=ChainFermi.cpp ChainPolymerized.cpp 
LADDER=LadderFermi.cpp LadderFree.cpp LadderFlux.cpp 
SQUARE=SquareFermi.cpp SquarePiFlux.cpp SquareACSL.cpp SquareFreeFlux.cpp SquareFreeHopping.cpp SquareJastrow.cpp 
TRIANGLE=TriangleFermi.cpp
HONEYCOMB=Honeycomb0pp.cpp
KAGOME=KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp

WF=$(CHAIN) $(LADDER) $(SQUARE) $(TRIANGLE) $(HONEYCOMB) $(KAGOME)
MC=MonteCarlo.cpp MCSystem.cpp Observable.cpp System.cpp IOSystem.cpp CreateSystem.cpp MCSim.cpp
MIN=VMCMinimization.cpp Interpolation.cpp VMCInterpolation.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp 
ANALYSE=Analyse.cpp AnalyseEnergy.cpp AnalyseChain.cpp AnalyseMagnetization.cpp AnalyseHoneycomb.cpp AnalyseLadder.cpp VMCMinimization.cpp Interpolation.cpp Directory.cpp 
IOFILES=Linux.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Gnuplot.cpp
OTHER=Lapack.cpp Parseur.cpp Fit.cpp

mc_SRCS=    mc.cpp    $(MC) $(IOFILES) $(OTHER) $(WF) 
min_SRCS=   min.cpp   $(MC) $(IOFILES) $(OTHER) $(WF) $(MIN)
check_SRCS= check.cpp $(MC) $(IOFILES) $(OTHER) $(WF) 
study_SRCS= study.cpp $(MC) $(IOFILES) $(OTHER) $(WF) $(MIN) $(ANALYSE)
load_SRCS=  load.cpp  $(MC) $(IOFILES) $(OTHER) $(WF) 

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
