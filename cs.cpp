/*!  @file cs.cpp */

#include "Chain.hpp"
#include "SquarePiFlux.hpp"
#include "SquareMu.hpp"
#include "SquareFermi.hpp"
#include "HoneycombSU4.hpp"

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
				HoneycombSU4 h(P);
				break;
			}
		case 4:
			{
				//double mu(0.0);
				//P.set("mu",mu);
				//if( std::abs(mu)<1e-14) { SquareFermi s(P); }
				//else { SquareMu s(P); }
				//SquareMu s(P);
				SquarePiFlux s(P);
				
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

