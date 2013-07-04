CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra 
#LDLIBS = -O3 -fopenmp -msse2

all:test

test.o:test.cpp Read.hpp Write.hpp Header.hpp Directory.hpp RST.hpp RSTfile.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp  Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Directory.o:Directory.cpp Directory.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

RSTfile.o:RSTfile.cpp RSTfile.hpp RST.hpp Write.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Read.o Write.o Header.o Directory.o RST.o RSTfile.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o test *.gch
