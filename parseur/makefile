CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra 

LDLIBS =  -pedantic -Wall -Wextra 

all:test 

test.o:test.cpp Parseur.hpp

test:test.o 

ref:
	doxygen Doxyfile

clean:
	rm *.o test
