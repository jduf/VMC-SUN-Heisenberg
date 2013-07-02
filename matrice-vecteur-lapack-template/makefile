CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test

test.o:test.cpp Matrix.hpp Lapack.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Lapack.o:Lapack.cpp Lapack.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Lapack.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch test
