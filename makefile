CXX = g++ 
CC = $(CXX) 

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2 
#CXXFLAGS = -pedantic -Wall -Wextra -larmadillo -DARMA_EXTRA_DEBUG
CXXFLAGS = -pedantic -Wall -Wextra -llapack -lblas -larmadillo
#CXXFLAGS = -O3 -fopenmp -msse2

#LIBFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS= -llapack -lblas -larmadillo
#LIBFLAGS = -O3 -fopenmp -msse2

all:test

test:test.o State.o System.o Chrono.o
	$(CXX) -o $@ $^ $(LDFLAGS)

State.o:State.cpp State.hpp System.hpp
	$(CXX) -c $(CXXFLAGS) $^

Chrono.o:Chrono.cpp Chrono.hpp
	$(CXX) -c $(CXXFLAGS) $^

test.o:test.cpp State.hpp
	$(CXX) -c $(CXXFLAGS) $^ 

clean:
	rm *.o test
