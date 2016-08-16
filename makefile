MACHINE=$(shell uname -n)_
ifneq (,$(filter $(MACHINE),ctmcpc33_ pink-floyd_))
	MACHINE=
	EXEC = check
	EXEC+= mc
	EXEC+= min
	#EXEC+= mcbi
	POSTPROCESS = cp $(EXEC) ../sim;
	EXEC+= study
	POSTPROCESS+= cp ../sim/study  ../

	CXX = g++ -std=c++14
	LAPACK = -llapack -lblas
	ERRORS = -Wall -Wextra -pedantic
	OPTION = -fopenmp
else
	EXEC= $(MACHINE)min
	EXEC+= $(MACHINE)mc
	POSTPROCESS= cp $(EXEC) ../sim/;
	EXEC+= $(MACHINE)study
	POSTPROCESS+= cp ../sim/$(MACHINE)study ../

	CXX = icpc -std=c++11
	LAPACK =# -Wl,--no-as-needed -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
	ERRORS = -Wall -Wextra -pedantic
	OPTION = -openmp -mkl -D MY_BLAS_ZDOTU

	#CXX = g++ -std=c++11
	#LAPACK = -Wl,--no-as-needed -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
	#ERRORS = -Wall -Wextra -pedantic
	#OPTION = -fopenmp
endif

CHAIN     = ChainFermi.cpp ChainPolymerized.cpp ChainFree.cpp ChainSAS.cpp
LADDER    = LadderFermi.cpp LadderFree.cpp LadderFreeFlux.cpp
SQUARE    = SquareFermi.cpp SquareFree.cpp SquareMu.cpp SquareDimerizedBar.cpp SquareT2x2.cpp SquareT3x2.cpp SquareT3x3.cpp SquareT4x2.cpp SquareT4x3.cpp SquareT4x4.cpp SquareLadder.cpp SquareVCS.cpp SquareFreeFlux.cpp SquarePiFlux.cpp SquareChiral.cpp SquareBox6.cpp SquareJastrow.cpp 
TRIANGLE  = TriangleFermi.cpp TriangleFree.cpp TrianglePlaquette.cpp TriangleAlternatingPlaquette.cpp TriangleMu.cpp TriangleT3x2.cpp TrianglePhi.cpp TriangleChiral.cpp
HONEYCOMB = HoneycombFermi.cpp HoneycombFree.cpp Honeycombp00.cpp Honeycomb0pp.cpp HoneycombPiFlux.cpp HoneycombChiral.cpp
KAGOME    = KagomeFermi.cpp KagomeFree.cpp KagomePlaquette3A.cpp KagomePlaquette3B.cpp KagomePlaquette6A.cpp  KagomePlaquette6B.cpp KagomeChiral.cpp KagomeChiralB.cpp KagomeVBC.cpp KagomePiHalfTriangle.cpp

WF         = $(CHAIN) $(LADDER) $(SQUARE) $(TRIANGLE) $(HONEYCOMB) $(KAGOME)
MONTECARLO = MonteCarlo.cpp MCSystem.cpp MCSim.cpp BiSystem.cpp
VMCMIN     = VMCMinimization.cpp Interpolation.cpp VMCInterpolation.cpp VMCExtract.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp VMCSystematic.cpp
ANALYSE    = Analyse.cpp AnalyseEnergy.cpp AnalyseMin.cpp AnalyseExtract.cpp AnalyseSystematic.cpp AnalyseChain.cpp AnalyseHoneycomb.cpp AnalyseLadder.cpp VMCMinimization.cpp VMCExtract.cpp VMCSystematic.cpp Interpolation.cpp Directory.cpp
IOFILES    = Linux.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Gnuplot.cpp
OTHER      = System.cpp IOSystem.cpp Observable.cpp CreateSystem.cpp Lapack.cpp Parseur.cpp Fit.cpp

$(MACHINE)check_SRCS = check.cpp $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
$(MACHINE)mc_SRCS    = mc.cpp    $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
$(MACHINE)mcbi_SRCS  = mcbi.cpp  $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
$(MACHINE)min_SRCS   = min.cpp   $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN)
$(MACHINE)study_SRCS = study.cpp $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN) $(ANALYSE)

#-----------------------------------------------------------------

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(ERRORS) $(OPTION) $(LAPACK) libcminpack.a

EXEC+=
$(EXEC)_SRCS=$(EXEC_SRCS)

BUILD=build

SRCS=$(wildcard *.cpp)

all: OPTION += -O3 -DNDEBUG
all:$(EXEC)
	$(POSTPROCESS)

debug: OPTION += -ggdb
debug:$(EXEC)
	$(POSTPROCESS)

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
