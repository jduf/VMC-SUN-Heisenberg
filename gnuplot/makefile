CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

all:test



test:test.o Gnuplot.o Write.o RST.o Header.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

Gnuplot.o:Gnuplot.cpp Gnuplot.hpp Write.hpp RST.hpp Header.hpp Time.hpp Matrix.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

test.o:test.cpp Gnuplot.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp  Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

	
clean:
	rm *.o *.gch 
