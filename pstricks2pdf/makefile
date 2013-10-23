CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3

all:pstricks2pdf

pstricks2pdf:pstricks2pdf.o Parseur.o Read.o Write.o RST.o Header.o PSTricks.o
	$(CXX) -o $@ $^ $(LDFLAGS)

pstricks2pdf.o:pstricks2pdf.cpp Directory.hpp Parseur.hpp PSTricks.hpp
	$(CXX) -c  $(CXXFLAGS) $^

PSTricks.o:PSTricks.cpp PSTricks.hpp Write.hpp Linux.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Parseur.o:Parseur.cpp Parseur.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

Read.o:Read.cpp Read.hpp  Header.hpp Matrix.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp  Header.hpp Matrix.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp  Matrix.hpp Vector.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp  
	$(CXX) -c  $(CXXFLAGS) $^


clean:
	rm *.o *.gch 
