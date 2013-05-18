CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDLIBS = -pedantic -Wall -Wextra
#LDLIBS = -O3 -fopenmp -msse2

all:test

test.o:test.cpp Chrono.hpp

Chrono.o:Chrono.cpp Chrono.hpp

test:test.o Chrono.o

ref:
	doxygen Doxyfile

clean:
	rm *.o
