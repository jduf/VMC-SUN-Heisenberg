EXEC=test

test_SRCS=test.cpp IOFiles.cpp RST.cpp Header.cpp Rand.cpp Gnuplot.cpp Fit.cpp

#-----------------------------------------------------------------

CXX = g++ -std=c++11

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
LAPACK = -llapack -lblas
OPTION = -O3 

BUILD=build

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

SRCS=$(wildcard *.cpp)

all:$(EXEC)

.SECONDEXPANSION:
$(EXEC): $$(patsubst %.cpp, $(BUILD)/%.o, $$($$@_SRCS)) 
	@echo Links $(notdir $^)
	$(CXX) -o $@ $^ libcminpack.a $(LDFLAGS) $(NOASSERT)

$(BUILD)/%.o:%.cpp
	@echo Creates $(notdir $@)
	$(CXX) -MD -c $(CXXFLAGS) $(NOASSERT)  $< -o $@

-include $(addprefix $(BUILD)/,$(SRCS:.cpp=.d))

clean:
	rm -f $(BUILD)/* test plot* data.dat
