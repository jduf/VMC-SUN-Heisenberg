/*!  @file min.cpp */

#include "Minimization.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	Minimization sim(P);
	sim.min(P.get<double>("xmax"));
}
