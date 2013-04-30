CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS = -pedantic -Wall -Wextra 
#LDLIBS = -O3 -fopenmp -msse2

all:test

test:test.o Rand.o Chrono.o

Rand.o:Rand.cpp Rand.hpp
Chono.o:Chrono.cpp Chrono.hpp

test.o:test.cpp Rand.hpp Chrono.hpp


clean:
	rm *.o test
