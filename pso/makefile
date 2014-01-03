CXX = g++
CC = $(CXX)

CXXFLAGS = -pedantic -Wall -Wextra -fopenmp

LDFLAGS =  -pedantic -Wall -Wextra -fopenmp

all:test test-openmp

test-openmp.o:test-openmp.cpp Time.hpp Write.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test.o:test.cpp Parseur.hpp Particle.hpp Swarm.hpp Rand.hpp Functions.hpp Time.hpp Write.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Rand.o:Rand.cpp Rand.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrix.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

Parseur.o:Parseur.cpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Parseur.o Rand.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

test-openmp:test-openmp.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch test
