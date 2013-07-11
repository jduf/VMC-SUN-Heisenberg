/*!  @file cs.cpp */

#include "Chain.hpp"
#include "Square.hpp"
#include "Honeycomb.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int lattice(P.get<unsigned int>("lattice"));
	switch(lattice){
		case 2:
			{
				Chain c(P);
				break;
			}
		case 3:
			{
				Honeycomb h(P);
				break;
			}
		case 4:
			{
				Square s(P);
				break;
			}
		default:
			{
				std::cerr<<"main : lattice need to be correctly defined"<<std::endl;
				std::cerr<<"       2 <-> chain"<<std::endl;
				std::cerr<<"       3 <-> honeycomb"<<std::endl;
				std::cerr<<"       4 <-> square"<<std::endl;
			}
	}
}

