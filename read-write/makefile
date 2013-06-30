CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra 
#LDLIBS = -O3 -fopenmp -msse2

all:test

test.o:test.cpp Read.hpp Write.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Array2D.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Array2D.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Array2D.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o test *.gch
