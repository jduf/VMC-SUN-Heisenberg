EXEC=test gnuplot2image

test_SRCS=test.cpp RST.cpp Gnuplot.cpp Write.cpp Header.cpp
gnuplot2image_SRCS=gnuplot2image.cpp Header.cpp RST.cpp Write.cpp Header.cpp Gnuplot.cpp Parseur.cpp

#-----------------------------------------------------------------

CXX = g++

NOASSERT = #-DNDEBUG
ERRORS = -Wall -Wextra -pedantic
OPTION = -O3 

BUILD=build

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(ERRORS) $(OPTION)

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
	rm -f $(BUILD)/* test bla* *.dat
