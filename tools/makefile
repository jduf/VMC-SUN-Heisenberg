CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3

all:diffdir

diffdir:diffdir.o Directory.o Parseur.o
	$(CXX) -o $@ $^ $(LDFLAGS)

diffdir.o:diffdir.cpp Directory.hpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Directory.o:Directory.cpp Directory.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

Parseur.o:Parseur.cpp Parseur.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

clean:
	rm *.o *.gch 
