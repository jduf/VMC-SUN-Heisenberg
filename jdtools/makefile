CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra 
#CXXFLAGS = -O3 -fopenmp -msse2

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra 
#LDLIBS = -O3 -fopenmp -msse2

all:jdwrite jdread jdfind

########
# jdfind
########
jdfind.o:jdfind.cpp Directory.hpp Read.hpp Write.hpp RST.hpp
	$(CXX) -c  $(CXXFLAGS) $^

jdfind:jdfind.o Header.o Read.o Write.o RST.o Directory.o
	$(CXX) -o $@ $^ $(LDFLAGS)

########
# jdread
########
jdread.o:jdread.cpp Header.hpp Read.hpp
	$(CXX) -c  $(CXXFLAGS) $^

jdread:jdread.o Header.o Read.o Write.o RST.o
	$(CXX) -o $@ $^ $(LDFLAGS)

#########
# jdwrite
#########
jdwrite.o:jdwrite.cpp Read.hpp Write.hpp Rewrite.hpp
	$(CXX) -c  $(CXXFLAGS) $^

jdwrite:jdwrite.o Read.o Write.o Header.o RST.o Rewrite.o
	$(CXX) -o $@ $^ $(LDFLAGS)

########
# commun
########
Read.o:Read.cpp Read.hpp Matrice.hpp Array2D.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Matrice.hpp Array2D.hpp Header.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrice.hpp Array2D.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Directory.o:Directory.cpp Directory.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp 
	$(CXX) -c  $(CXXFLAGS) $^



ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch 
