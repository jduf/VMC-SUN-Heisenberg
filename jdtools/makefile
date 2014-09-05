EXEC=jdinfo jdread jdtags jdfind

jdinfo_SRCS = jdinfo.cpp Directory.cpp Header.cpp Parseur.cpp IOFiles.cpp RST.cpp RSTFile.cpp
jdread_SRCS = jdread.cpp Directory.cpp Header.cpp Parseur.cpp IOFiles.cpp RST.cpp 
jdtags_SRCS = jdtags.cpp Directory.cpp Header.cpp Parseur.cpp IOFiles.cpp RST.cpp RSTFile.cpp
jdfind_SRCS = jdfind.cpp Directory.cpp Header.cpp Parseur.cpp IOFiles.cpp RST.cpp RSTFile.cpp
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
	rm -f $(BUILD)/*
