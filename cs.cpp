/*!  @file cs.cpp */

#include "CreateSystem.hpp"
#include "Parseur.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem CS(P);
	CS.create();
	CS.save();
}
