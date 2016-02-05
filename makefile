MACHINE=$(shell uname -n)_
ifneq (,$(filter $(MACHINE),ctmcpc33_ pink-floyd_))
	MACHINE=
	EXEC = min
	EXEC+= mc
	#EXEC+= mcbi
	#EXEC+= load
	EXEC+= check
	POSTPROCESS= cp $(EXEC) ../sim;
	#EXEC+= study
	#POSTPROCESS+= cp study ..

	CXX = g++ -std=c++11
	LAPACK = -llapack -lblas
	ERRORS = -Wall -Wextra -pedantic
	OPTION = -fopenmp
else
	EXEC = $(MACHINE)min
	POSTPROCESS= cp $(EXEC) ../sim;
	EXEC+= $(MACHINE)study
	POSTPROCESS+= cp $(MACHINE)study ..

	CXX = icpc -std=c++11
	LAPACK =# -Wl,--no-as-needed -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
	ERRORS = -Wall -Wextra -pedantic
	OPTION = -openmp -mkl -D MY_BLAS_ZDOTU

	#CXX = g++ -std=c++11
	#LAPACK = -Wl,--no-as-needed -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
	#ERRORS = -Wall -Wextra -pedantic
	#OPTION = -fopenmp
endif

CHAIN     = ChainFermi.cpp ChainPolymerized.cpp ChainFree.cpp
LADDER    = LadderFermi.cpp LadderFree.cpp LadderFreeFlux.cpp
SQUARE    = SquareFermi.cpp SquareFree.cpp SquareMu.cpp SquareFreeFlux.cpp SquarePiFlux.cpp SquareChiral.cpp SquareJastrow.cpp
TRIANGLE  = TriangleFermi.cpp TriangleFree.cpp TrianglePlaquette.cpp TriangleMu.cpp TrianglePhi.cpp TriangleChiral.cpp
HONEYCOMB = HoneycombFermi.cpp HoneycombFree.cpp HoneycombPlaquette.cpp Honeycomb0pp.cpp HoneycombPiFlux.cpp HoneycombChiral.cpp
KAGOME    = KagomeFermi.cpp KagomeDirac.cpp KagomeVBC.cpp

WF         = $(CHAIN) $(LADDER) $(SQUARE) $(TRIANGLE) $(HONEYCOMB) $(KAGOME)
MONTECARLO = MonteCarlo.cpp MCSystem.cpp Observable.cpp System.cpp IOSystem.cpp CreateSystem.cpp MCSim.cpp BiSystem.cpp
VMCMIN     = VMCMinimization.cpp Interpolation.cpp VMCInterpolation.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp VMCSystematic.cpp
ANALYSE    = Analyse.cpp AnalyseEnergy.cpp AnalyseChain.cpp AnalyseHoneycomb.cpp AnalyseMagnetization.cpp AnalyseMin.cpp VMCMinimization.cpp Interpolation.cpp Directory.cpp
IOFILES    = Linux.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Gnuplot.cpp
OTHER      = Lapack.cpp Parseur.cpp Fit.cpp

$(MACHINE)mc_SRCS    = mc.cpp    $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
$(MACHINE)mcbi_SRCS  = mcbi.cpp  $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
$(MACHINE)min_SRCS   = min.cpp   $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN)
$(MACHINE)check_SRCS = check.cpp $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
$(MACHINE)study_SRCS = study.cpp $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN) $(ANALYSE)
$(MACHINE)load_SRCS  = load.cpp  $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)

#-----------------------------------------------------------------

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(ERRORS) $(OPTION) $(LAPACK) libcminpack.a

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
