CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS = -pedantic -Wall -Wextra
#LDLIBS = -O3 -fopenmp -msse2

all:test

Complex.o:Complex.cpp Complex.hpp
test.o:test.cpp

test:test.o Complex.o 


clear:
	rm *.o
