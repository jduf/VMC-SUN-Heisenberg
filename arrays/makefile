CXX = g++
CC = $(CXX)

CXXFLAGS = -pedantic -Wall -Wextra -llapack -lblas

LDLIBS = -pedantic -Wall -Wextra -llapack -lblas

all:test 

test.o:test.cpp Array1D.hpp Array2D.hpp Read.hpp Write.hpp Matrix.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp Matrix.hpp Array2D.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Matrix.hpp Array2D.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrix.hpp Array2D.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Read.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o test *.gch
