CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

all:run

run:run.o Directory.o Parseur.o RSTfile.o Write.o Header.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

run.o:run.cpp Linux.hpp Parseur.hpp RSTfile.hpp RST.hpp Write.hpp Header.hpp Matrix.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Directory.o:Directory.cpp Directory.hpp 
	$(CXX) -c  $(CXXFLAGS) $^
	
Parseur.o:Parseur.cpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RSTfile.o:RSTfile.cpp RSTfile.hpp RST.hpp Write.hpp
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp RST.hpp Write.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp  Header.hpp Matrix.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrix.hpp Time.hpp
	$(CXX) -c  $(CXXFLAGS) $^

clean:
	rm *.o *.gch 
