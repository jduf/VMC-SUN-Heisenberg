#include "CreateSystem.hpp"

CreateSystem::CreateSystem(Parseur& P):
	param_(false),
	ref_(3,0),
	Sr_(NULL),
	Sc_(NULL),
	valid_(true)
{
	std::string wf(P.get<std::string>("wf"));
	param_.set("N",P.get<unsigned int>("N"));
	param_.set("n",P.get<unsigned int>("n"));
	param_.set("m",P.get<unsigned int>("m"));

	
	if( wf == "chainfermi" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "chaindimerized" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 1;
		param_.set("delta",P.get<unsigned int>("delta"));
	}
	if( wf == "chaintrimerized" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 0;
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
		param_.set("mu",P.get<double>("mu"));
	}
	if( wf == "trianglephi" ){
		ref_(0) = 3;
		ref_(1) = 2;
		ref_(2) = 2;
		param_.set("phi",P.get<double>("phi"));
	}
	if( wf == "trianglejastrow" ){
		ref_(0) = 3;
		ref_(1) = 0;
	}

	if( wf == "squaref_ermi" ){
		ref_(0) = 4;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "squaremu" ){
		ref_(0) = 4;
		ref_(1) = 1;
		ref_(2) = 1;
		param_.set("mu",P.get<double>("mu"));
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
		param_.set("phi",P.get<double>("phi"));
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

	param_.set("ref",ref_);

	if(P.status()){
		valid_ = false; 
		std::cout<<"CreateSystem::CreateSystem(P) : some parameters are missing, can't create the system"<<std::endl;
	}
}

CreateSystem::~CreateSystem(){
	if(Sr_){delete Sr_;}
	if(Sc_){delete Sc_;}
}

void CreateSystem::create(){
	if(valid_){
		std::cout<<"bla"<<std::endl;
		switch(ref_(0)){
			case 2:
				{
					switch(ref_(1)){
						case 1:
							{
								switch(ref_(2)) {
									case 0:{Sr_ = new ChainFermi(param_);}break;
									case 1:{Sr_ = new ChainDimerized(param_);}break;
									case 3:{Sr_ = new ChainTrimerized(param_);}break;
									default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								}
							}break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			case 3:
				{
					switch(ref_(1)){
						case 0:{Sr_ = new TriangleJastrow(param_);}break;
						case 1:
							   {
								   switch(ref_(2)){
									   case 0:{Sr_ = new TriangleFermi(param_);}break;
									   case 1:{Sr_ = new TriangleMu(param_);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
								   }
							   }break;
						case 2:
							   {
								   switch(ref_(2)){
									   case 4:{Sc_ = new TrianglePhi(param_);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
								   }
							   }break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			case 4:
				{
					switch(ref_(1)){
						case 0:{Sr_ = new SquareJastrow(param_);}break;
						case 1:
							   {
								   switch(ref_(2)){
									   case 0:{Sr_ = new SquareFermi(param_);}break;
									   case 1:{Sr_ = new SquareMu(param_);}break;
									   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								   }
							   }break;
						case 2:
							   {
								   switch(ref_(2)){
									   case 3:{Sc_ = new SquarePiFlux(param_);}break;
									   case 4:{Sc_ = new SquareSU2PhiFlux(param_);}break;
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
									case 0:{Sr_ = new HoneycombSU3(param_);}break;
									case 1:{Sr_ = new HoneycombSU4(param_);}break;
									default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
								}
							}break;
						default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
					}
				}break;
			default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
		}
	}
}

void CreateSystem::save(){
	if(Sr_){Sr_->save();}
	if(Sc_){Sc_->save();}
}
