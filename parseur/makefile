CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra 

LDLIBS =  -pedantic -Wall -Wextra 

all:test 

test:test.o Parseur.o

test.o:test.cpp Parseur.hpp

Parseur.o:Parseur.cpp Parseur.hpp

ref:
	doxygen Doxyfile

clean:
	rm *.o test
