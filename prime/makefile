CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS =  -pedantic -Wall -Wextra -llapack -lblas 
#LDLIBS = -O3 -fopenmp -msse2

all:test


test.o:test.cpp Prime.hpp

test:test.o 


clean:
	rm *.o test
