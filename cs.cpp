/*!  @file cs.cpp */

#include "ChainFermi.hpp"
#include "ChainDimerized.hpp"

#include "SquareFermi.hpp"
#include "SquareMu.hpp"
#include "SquarePiFlux.hpp"
#include "SquareSU2PhiFlux.hpp"
#include "SquareJastrow.hpp"

#include "HoneycombSU4.hpp"
#include "HoneycombSU3.hpp"

#include "TriangleFermi.hpp"
#include "TriangleMu.hpp"
#include "TrianglePhi.hpp"
#include "TriangleJastrow.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string wf(P.get<std::string>("wf"));
	if( wf == "chainfermi" ){
		ChainFermi c(P);
		c.save();
	}
	if( wf == "chaindimerized" ){
		ChainDimerized c(P);
		c.save();
	}

	if( wf == "squarefermi" ){
		SquareFermi s(P);
		s.save();
	}
	if( wf == "squaremu" ){
		SquareMu s(P);
		s.save();
	}
	if( wf == "squarecsl" ){
		SquarePiFlux s(P);
		s.save();
	}
	if( wf == "squarephi" ){
		SquareSU2PhiFlux s(P);
		s.save();
	}
	if( wf == "squarejastrow" ){
		SquareJastrow s(P);
		s.save();
	}

	if( wf == "trianglefermi" ){
		TriangleFermi s(P);
		s.save();
	}
	if( wf == "trianglemu" ){
		TriangleMu s(P);
		s.save();
	}
	if( wf == "trianglephi" ){
		TrianglePhi s(P);
		s.save();
	}
	if( wf == "trianglejastrow" ){
		TriangleJastrow s(P);
		s.save();
	}

	if( wf == "honeycombsu4" ){
		HoneycombSU4 h(P);
		h.save();
	}
	if( wf == "honeycombsu3" ){
		HoneycombSU3 h(P);
		h.save();
	}
}

