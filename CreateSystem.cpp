#include "CreateSystem.hpp"

CreateSystem::CreateSystem(Parseur& P):
	param_(false),
	Sr_(NULL),
	Sc_(NULL),
	ref_(3,0)
{
	parse(P,param_);
	init(param_);
}

CreateSystem::CreateSystem(Container const& C):
	param_(false),
	Sr_(NULL),
	Sc_(NULL),
	ref_(C.get<Vector<unsigned int> >("ref"))
{
	init(C);
}

CreateSystem::~CreateSystem(){
	if(Sr_){delete Sr_;}
	if(Sc_){delete Sc_;}
}

void CreateSystem::parse(Parseur& P, Container& C) {
	std::string wf(P.get<std::string>("wf"));
	C.set("N",P.get<unsigned int>("N"));
	C.set("n",P.get<unsigned int>("n"));
	C.set("m",P.get<unsigned int>("m"));

	if( wf == "chainfermi" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "chaindimerized" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 1;
	}
	if( wf == "chainpolymerized" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 2;
	}

	if( wf == "trianglefermi" ){
		ref_(0) = 3;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "trianglemu" ){
		ref_(0) = 3;
		ref_(1) = 1;
		ref_(2) = 1;
	}
	if( wf == "trianglephi" ){
		ref_(0) = 3;
		ref_(1) = 2;
		ref_(2) = 2;
	}
	if( wf == "trianglejastrow" ){
		ref_(0) = 3;
		ref_(1) = 0;
	}

	if( wf == "squarefermi" ){
		ref_(0) = 4;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "squaremu" ){
		ref_(0) = 4;
		ref_(1) = 1;
		ref_(2) = 1;
	}
	if( wf == "squarecsl" ){
		ref_(0) = 4;
		ref_(1) = 1;
		ref_(2) = 2;
	}
	if( wf == "squarephi" ){
		ref_(0) = 4;
		ref_(1) = 2;
		ref_(2) = 3;
	}
	if( wf == "squarejastrow" ){
		ref_(0) = 4;
		ref_(1) = 0;
		ref_(2) = 4;
	}

	if( wf == "honeycombsu3" ){
		ref_(0) = 6;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "honeycombsu4" ){
		ref_(0) = 6;
		ref_(1) = 1;
		ref_(2) = 0;
	}

	if(ref_(0)==2){
		C.set("a",P.get<unsigned int>("a"));
	}

	C.set("ref",ref_);
}

void CreateSystem::init(Container const& C){
		switch(ref_(0)){
			case 2:
				{
					switch(ref_(1)){
						case 1:
							{
								switch(ref_(2)){
									case 0:{Sr_ = new ChainFermi(C);}break;
									case 1:{Sr_ = new ChainDimerized(C);}break;
									case 2:{Sr_ = new ChainPolymerized(C);}break;
									default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								}
							}break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			case 3:
				{
					switch(ref_(1)){
						case 0:{Sr_ = new TriangleJastrow(C);}break;
						case 1:
							   {
								   switch(ref_(2)){
									   case 0:{Sr_ = new TriangleFermi(C);}break;
									   case 1:{Sr_ = new TriangleMu(C);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
								   }
							   }break;
						case 2:
							   {
								   switch(ref_(2)){
									   case 4:{Sc_ = new TrianglePhi(C);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
								   }
							   }break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			case 4:
				{
					switch(ref_(1)){
						case 0:{Sr_ = new SquareJastrow(C);}break;
						case 1:
							   {
								   switch(ref_(2)){
									   case 0:{Sr_ = new SquareFermi(C);}break;
									   case 1:{Sr_ = new SquareMu(C);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								   }
							   }break;
						case 2:
							   {
								   switch(ref_(2)){
									   case 3:{Sc_ = new SquarePiFlux(C);}break;
									   case 4:{Sc_ = new SquareSU2PhiFlux(C);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								   }
							   }break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			case 6:
				{
					switch(ref_(1)){
						case 1:
							{
								switch(ref_(2)){
									case 0:{Sr_ = new HoneycombSU3(C);}break;
									case 1:{Sr_ = new HoneycombSU4(C);}break;
									default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								}
							}break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
		}
}

void CreateSystem::save(){
	if(Sr_){Sr_->save();}
	if(Sc_){Sc_->save();}
}

void CreateSystem::get_param(Container& param){
	if(Sr_) {Sr_->get_input(param);}
	if(Sc_) {Sc_->get_input(param);}
}

void CreateSystem::get_input(Container& input){
	if(Sr_) {Sr_->get_input(input);}
	if(Sc_) {Sc_->get_input(input);}
}

bool CreateSystem::use_complex(){
	if(ref_(1) == 1){ return false; }
	else { return true; }
}
