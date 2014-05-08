EXEC=test

test_SRCS=test.cpp Rand.cpp Write.cpp Gnuplot.cpp RST.cpp Header.cpp 

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
