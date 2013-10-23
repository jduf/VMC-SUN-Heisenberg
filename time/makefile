CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

all:test

test:test.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

test.o:test.cpp Time.hpp 
	$(CXX) -c  $(CXXFLAGS) $^
	
clean:
	rm *.o *.gch 
