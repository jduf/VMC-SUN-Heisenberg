CXX = g++
CC = $(CXX)

CXXFLAGS = -pedantic -Wall -Wextra 

LDFLAGS = -pedantic -Wall -Wextra 

all:test

test.o:test.cpp Read.hpp Write.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp  Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o test *.gch
