EXEC=test

test_SRCS=test.cpp Fit.cpp Lapack.cpp Read.cpp Header.cpp RST.cpp Write.cpp Gnuplot.cpp

#-----------------------------------------------------------------

CXX = g++

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
	$(CXX) -o $@ $^ liblevmar.a $(LDFLAGS) $(NOASSERT)

$(BUILD)/%.o:%.cpp
	@echo Creates $(notdir $@)
	$(CXX) -MD -c $(CXXFLAGS) $(NOASSERT)  $< -o $@

-include $(addprefix $(BUILD)/,$(SRCS:.cpp=.d))

clean:
	rm -f $(BUILD)/* test fit*
