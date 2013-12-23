CXX = g++
CC = $(CXX)

CXXFLAGS = -pedantic -Wall -Wextra 

LDFLAGS = -pedantic -Wall -Wextra 

all:evalnote

evalnote.o:evalnote.cpp Read.hpp Write.hpp Header.hpp Matrix.hpp Vector.hpp Parseur.hpp Gnuplot.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Gnuplot.o:Gnuplot.cpp Gnuplot.hpp Write.hpp RST.hpp Header.hpp Time.hpp Matrix.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp  Header.hpp Matrix.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp  Header.hpp Matrix.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp  Matrix.hpp Vector.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^

Parseur.o:Parseur.cpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

evalnote:evalnote.o Read.o Write.o Header.o RST.o Parseur.o Gnuplot.o
	$(CXX) -o $@ $^ $(LDFLAGS)

ref:
	doxygen Doxyfile

clean:
	rm *.o test *.gch
