/*! @file min.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);

	IOFiles read(P.get<std::string>("sim"),false);

	Vector<double> tmp;
	read>>tmp;
	System s(read);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);

	cs.lattice("./","lattice");
}
