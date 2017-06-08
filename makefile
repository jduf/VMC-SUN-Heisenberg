EXEC=
EXEC+= mc
EXEC+= min
EXEC+= check
EXEC+= study
#EXEC+= mcbi

CXX = g++ -std=c++14
LAPACK = -llapack -lblas
OPTION = -fopenmp
ERRORS = -Wall -Wextra -pedantic

#CXX = icpc -std=c++14
#LAPACK =# -Wl,--no-as-needed -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
#OPTION = -qopenmp -mkl -D MY_BLAS_ZDOTU -no-inline-max-size -no-inline-max-total-size
ERRORS = -Wall -Wextra -pedantic

JDLIB = lib/lib
BUILD = build
#ADAPT DEPENDENCIES -----------------------------------------------------------------

CHAIN     = ChainFermi.cpp ChainPolymerized.cpp ChainFree.cpp ChainSAS.cpp
LADDER    = LadderFermi.cpp LadderFree.cpp LadderDimerA.cpp LadderDimerB.cpp LadderSquarePlaquetteA.cpp LadderSquarePlaquetteB.cpp LadderSquarePlaquetteC.cpp LadderRectangularPlaquetteA.cpp LadderRectangularPlaquetteB.cpp LadderRectangularPlaquetteC.cpp LadderRectangularPlaquetteD.cpp
SQUARE    = SquareFermi.cpp SquareFree.cpp SquareMu.cpp SquareDimerizedBar.cpp SquareT2x2.cpp SquareT3x2.cpp SquareT3x3.cpp SquareT4x2.cpp SquareT4x3.cpp SquareT4x4.cpp SquareLadder.cpp SquareVCS.cpp SquareMuT2x1.cpp SquareMuT2x2.cpp SquareFacingDimers.cpp SquareAlternatingDimers.cpp Squarek2MuPi.cpp SquareFreeFlux.cpp SquarePiFlux.cpp SquareChiral.cpp SquareBox6.cpp SquareJastrow.cpp
TRIANGLE  = TriangleFermi.cpp TriangleFree.cpp TrianglePlaquette.cpp TriangleAlternatingPlaquette.cpp TriangleMu.cpp TriangleT3x2.cpp TrianglePhi.cpp TriangleChiral.cpp TriangleChiralSG.cpp
HONEYCOMB = HoneycombFermi.cpp HoneycombFree.cpp Honeycombp00.cpp Honeycomb0pp.cpp HoneycombPiFlux.cpp HoneycombChiral.cpp
KAGOME    = KagomeFermi.cpp KagomeFree.cpp KagomePlaquette3A.cpp KagomePlaquette3B.cpp KagomePlaquette6A.cpp  KagomePlaquette6B.cpp KagomeChiral.cpp KagomeChiralB.cpp KagomeVBC.cpp KagomePiHalfTriangle.cpp

WF         = $(CHAIN) $(LADDER) $(SQUARE) $(TRIANGLE) $(HONEYCOMB) $(KAGOME)
MONTECARLO = MonteCarlo.cpp MCSystem.cpp MCSim.cpp BiSystem.cpp
VMCMIN     = VMCMinimization.cpp Interpolation.cpp VMCInterpolation.cpp VMCExtract.cpp PSO.cpp VMCPSO.cpp MCParticle.cpp VMCSystematic.cpp ACiD.cpp VMCACiD.cpp VMCAnalyse.cpp Directory.cpp
ANALYSE    = Analyse.cpp AnalyseEnergy.cpp AnalyseMin.cpp AnalyseExtract.cpp AnalyseSystematic.cpp AnalyseSystematicCoupling.cpp AnalyseChain.cpp AnalyseHoneycomb.cpp AnalyseLadder.cpp Interpolation.cpp
IOFILES    = Linux.cpp IOFiles.cpp Header.cpp RST.cpp RSTFile.cpp PSTricks.cpp Gnuplot.cpp
OTHER      = System.cpp IOSystem.cpp Observable.cpp CreateSystem.cpp Lapack.cpp Parseur.cpp Fit.cpp

mc_SRCS    = mc.cpp    $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)
min_SRCS   = min.cpp   $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN)
check_SRCS = check.cpp $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN) $(ANALYSE)
study_SRCS = study.cpp $(MONTECARLO) $(IOFILES) $(OTHER) $(WF) $(VMCMIN) $(ANALYSE)
mcbi_SRCS  = mcbi.cpp  $(MONTECARLO) $(IOFILES) $(OTHER) $(WF)

#DO NOT CHANGE BELOW ----------------------------------------------------------------

VPATH = $(JDLIB) .
MYLIB = $(addprefix -I,$(VPATH))
-include $(JDLIB)/config.mk

CXXFLAGS = $(MYLIB) $(ERRORS) $(OPTION) $(CONFIG)
LDFLAGS  = $(MYLIB) $(ERRORS) $(OPTION) $(LAPACK) $(JDLIB)/libcminpack.a

all: OPTION += -O3 -DNDEBUG
all:$(EXEC)

debug: OPTION += -ggdb
debug:$(EXEC)

.SECONDEXPANSION:
$(EXEC): $$(patsubst %.cpp, $(BUILD)/%.o, $$($$@_SRCS))
	@echo Links $(notdir $^) 
	$(CXX) $^ $(LDFLAGS) -o $@

$(BUILD)/%.o:%.cpp
	@echo Creates $(notdir $@)
	$(CXX) $< $(CXXFLAGS) -MMD -c -o $@

-include $(wildcard $(BUILD)/*.d)

clean:
	rm -f $(BUILD)/* $(EXEC)

ref:
	@echo Create the documentation
	@doxygen doxygen/Doxyfile > doxygen/log
	@firefox doxygen/html/files.html &
