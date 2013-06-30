CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS =  -pedantic -Wall -Wextra -llapack -lblas 
#LDLIBS = -O3 -fopenmp -msse2

all:test


Prime.o:Prime.cpp Prime.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Prime.o
	$(CXX) -o $@ $^ $(LDFLAGS)


clean:
	rm *.o test
