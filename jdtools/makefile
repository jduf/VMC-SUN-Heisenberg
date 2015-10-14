EXEC=jdread jdhtml

jdread_SRCS = jdread.cpp Header.cpp IOFiles.cpp RST.cpp Linux.cpp
jdhtml_SRCS = jdhtml.cpp Header.cpp IOFiles.cpp RST.cpp RSTFile.cpp Linux.cpp 

#-----------------------------------------------------------------

CXX = g++ -std=c++11

ERRORS = -Wall -Wextra -pedantic
OPTION = 

BUILD=build

CXXFLAGS = $(ERRORS) $(OPTION)
LDFLAGS  = $(ERRORS) $(OPTION)

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
	rm -f $(BUILD)/*
