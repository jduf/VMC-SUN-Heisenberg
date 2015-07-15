EXEC=test

test_SRCS=test.cpp RST.cpp Header.cpp IOFiles.cpp Interpolation.cpp Gnuplot.cpp Lapack.cpp

#-----------------------------------------------------------------

CXX = g++ -std=c++11

ERRORS = -Wall -Wextra -pedantic
OPTION = -fopenmp
LAPACK = -llapack -lblas

BUILD=build

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(LAPACK) $(ERRORS) $(OPTION)

SRCS=$(wildcard *.cpp)

all: OPTION += -O3 -DNDEBUG
all:$(EXEC)

debug: OPTION += -ggdb 
debug:$(EXEC)

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
