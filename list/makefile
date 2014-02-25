CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test

test.o:test.cpp List.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch test
