/*!  @file cs.cpp */

#include "Chain.hpp"

#include "SquareFermi.hpp"
#include "SquareMu.hpp"
#include "SquarePiFlux.hpp"
#include "SquareSU2PhiFlux.hpp"
#include "SquareSU2AF.hpp"

#include "HoneycombSU4.hpp"

#include "TriangleFermi.hpp"
#include "TriangleMu.hpp"
#include "TrianglePhi.hpp"

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
	if( wf == "phi" ){
		SquareSU2PhiFlux s(P);
	}
	if( wf == "AF" ){
		SquareSU2AF s(P);
	}
	if(wf == "triangle" ){
		TriangleFermi s(P);
	}
	if(wf == "trianglemu" ){
		TriangleMu s(P);
	}
	if(wf == "trianglephi" ){
		TrianglePhi s(P);
	}
}

