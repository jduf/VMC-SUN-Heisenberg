CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3

all:test

test:test.o
	$(CXX) -o $@ $^ $(LDFLAGS)

test.o:test.cpp Linux.hpp
	$(CXX) -c  $(CXXFLAGS) $^

clean:
	rm *.o *.gch 
