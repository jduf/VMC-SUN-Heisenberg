EXEC=test gnuplot2image

test_SRCS=test.cpp Gnuplot.cpp IOFiles.cpp Header.cpp RST.cpp 
#gnuplot2image_SRCS=gnuplot2image.cpp Gnuplot.cpp Parseur.cpp IOFiles.cpp Header.cpp RST.cpp 
gnuplot2image_SRCS=gnuplot2image.cpp 

#-----------------------------------------------------------------

CXX = g++ -std=c++11

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
