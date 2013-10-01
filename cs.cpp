/*!  @file cs.cpp */

#include "Chain.hpp"
#include "SquarePiFlux.hpp"
#include "SquareMu.hpp"
#include "SquareFermi.hpp"
#include "HoneycombSU4.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string wf(P.get<std::string>("wf"));
	if( wf == "chain" ){
		Chain c(P);
	}
	if( wf == "fermi" ){
		SquareFermi s(P);
	}
	if( wf == "mu" ){
		SquareMu s(P);
	}
	if( wf == "csl" ){
		SquarePiFlux s(P);
	}
	if( wf == "honeycomb" ){
		HoneycombSU4 h(P);
	}
}

