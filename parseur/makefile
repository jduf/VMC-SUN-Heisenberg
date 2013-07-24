CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra 

LDLIBS =  -pedantic -Wall -Wextra 

all:test 

test:test.o Parseur.o
	$(CXX) -o $@ $^ $(LDFLAGS)

test.o:test.cpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Parseur.o:Parseur.cpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

ref:
	doxygen Doxyfile

clean:
	rm *.o test
