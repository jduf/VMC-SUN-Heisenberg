CXX = g++ -std=c++11
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test

test.o:test.cpp Matrix.hpp Lapack.hpp Read.hpp Gnuplot.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Lapack.o:Lapack.cpp Lapack.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrix.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

Gnuplot.o:Gnuplot.cpp Gnuplot.hpp Write.hpp RST.hpp Header.hpp Time.hpp Matrix.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Lapack.o Read.o Write.o Header.o RST.o Gnuplot.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch test
