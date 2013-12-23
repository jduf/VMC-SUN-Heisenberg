CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS = -pedantic -Wall -Wextra 
#LDLIBS = -O3 -fopenmp -msse2

all:test

test:test.o Rand.o 

Rand.o:Rand.cpp Rand.hpp

test.o:test.cpp Rand.hpp Time.hpp


clean:
	rm *.o test
