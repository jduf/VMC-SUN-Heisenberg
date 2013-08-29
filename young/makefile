CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp

all:young

young:young.o YoungTableau.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

young.o:young.cpp YoungTableau.hpp 
	$(CXX) -c  $(CXXFLAGS) $^
	
YoungTableau.o:YoungTableau.cpp YoungTableau.hpp 
	$(CXX) -c  $(CXXFLAGS) $^
	
clean:
	rm *.o *.gch 
