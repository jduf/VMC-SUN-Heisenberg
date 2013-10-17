CXX = g++
CC = $(CXX)

#CXXFLAGS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
CXXFLAGS = -pedantic -Wall -Wextra -O3

#LDLIBS = -pedantic -Wall -Wextra -O3 -fopenmp -msse2
LDFLAGS = -pedantic -Wall -Wextra -O3

all:jdwrite jdread jdfind jdtags jdinfo

########
# jdfind
########
jdfind:jdfind.o Header.o Read.o Write.o RST.o RSTfile.o Directory.o Parseur.o
	$(CXX) -o $@ $^ $(LDFLAGS)

jdfind.o:jdfind.cpp Directory.hpp Read.hpp Write.hpp RSTfile.hpp Parseur.hpp Linux.hpp
	$(CXX) -c  $(CXXFLAGS) $^

########
# jdread
########
jdread:jdread.o Header.o Read.o Write.o RST.o Parseur.o
	$(CXX) -o $@ $^ $(LDFLAGS)

jdread.o:jdread.cpp Header.hpp Read.hpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

#########
# jdwrite
#########
jdwrite:jdwrite.o Read.o Write.o Header.o RST.o Rewrite.o
	$(CXX) -o $@ $^ $(LDFLAGS)

jdwrite.o:jdwrite.cpp Read.hpp Write.hpp Rewrite.hpp Linux.hpp
	$(CXX) -c  $(CXXFLAGS) $^

#########
# jdtags
#########
jdtags:jdtags.o Read.o Write.o Header.o RST.o RSTfile.o Parseur.o
	$(CXX) -o $@ $^ $(LDFLAGS)

jdtags.o:jdtags.cpp Read.hpp Write.hpp RSTfile.hpp Linux.hpp Matrix.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

##########
# jdinfo
##########
jdinfo:jdinfo.o Header.o Read.o Write.o RST.o RSTfile.o Directory.o Parseur.o
	$(CXX) -o $@ $^ $(LDFLAGS)

jdinfo.o:jdinfo.cpp Directory.hpp Read.hpp Write.hpp RSTfile.hpp Parseur.hpp
	$(CXX) -c  $(CXXFLAGS) $^

########
# commun
########
Read.o:Read.cpp Read.hpp Matrix.hpp Header.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Write.o:Write.cpp Write.hpp Matrix.hpp Header.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Header.o:Header.cpp Header.hpp RST.hpp Matrix.hpp Time.hpp Vector.hpp
	$(CXX) -c  $(CXXFLAGS) $^

Directory.o:Directory.cpp Directory.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

Parseur.o:Parseur.cpp Parseur.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

RST.o:RST.cpp RST.hpp 
	$(CXX) -c  $(CXXFLAGS) $^

RSTfile.o:RSTfile.cpp RSTfile.hpp Write.hpp
	$(CXX) -c  $(CXXFLAGS) $^

ref:
	doxygen Doxyfile

clean:
	rm *.o *.gch 
