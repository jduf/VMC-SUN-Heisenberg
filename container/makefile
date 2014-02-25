
CXX = g++
CC = $(CXX)

CXXFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

LDFLAGS =  -pedantic -Wall -Wextra -llapack -lblas 

all:test

Read.o:Read.cpp Read.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^
	
Parseur.o:Parseur.cpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Container.o:Container.cpp Container.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrix.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

test.o:test.cpp Container.hpp Read.hpp Header.hpp RST.hpp Parseur.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

test:test.o Read.o Header.o RST.o Parseur.o Container.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch test
