EXEC=jdinfo jdread jdtags jdwrite jdfind

jdinfo_SRCS= jdinfo.cpp Directory.cpp Header.cpp Parseur.cpp Read.cpp RST.cpp Write.cpp RSTFile.cpp
jdread_SRCS= jdread.cpp Directory.cpp Header.cpp Parseur.cpp Read.cpp RST.cpp Write.cpp
jdtags_SRCS= jdtags.cpp Directory.cpp Header.cpp Parseur.cpp Read.cpp RST.cpp Write.cpp RSTFile.cpp
jdfind_SRCS= jdfind.cpp Directory.cpp Header.cpp Parseur.cpp Read.cpp RST.cpp Write.cpp RSTFile.cpp
jdwrite_SRCS=jdwrite.cpp Directory.cpp Header.cpp Parseur.cpp Read.cpp Rewrite.cpp RST.cpp Write.cpp
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
	rm -f $(BUILD)/*
