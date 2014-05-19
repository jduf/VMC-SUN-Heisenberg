#include "CreateSystem.hpp"

CreateSystem::CreateSystem(Parseur& P):
	N_(P.get<unsigned int>("N")),
	n_(P.get<unsigned int>("n")),
	m_(P.get<unsigned int>("m")),
	param_(0),
	bc_(P.get<int>("bc")),
	ready_(false),
	ref_(3,0),
	RGL_(NULL),
	CGL_(NULL)
{
	parse(P);
	create();
}

CreateSystem::CreateSystem(CreateSystem const& cs, double param):
	N_(cs.N_),
	n_(cs.n_),
	m_(cs.m_),
	param_(param),
	bc_(cs.bc_),
	ready_(cs.ready_),
	ref_(cs.ref_),
	RGL_(NULL),
	CGL_(NULL)
{
	create();
	if(RGL_){ready_ = RGL_->create(param);}
	if(CGL_){ready_ = CGL_->create(param);}
}

CreateSystem::~CreateSystem(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
}

template<>
Matrix<double> CreateSystem::get_EVec() const { 
	if(RGL_){return RGL_->get_EVec();}
	return 0;
}

template<>
Matrix<std::complex<double> > CreateSystem::get_EVec() const { 
	if(CGL_){return CGL_->get_EVec();}
	return 0;
}

void CreateSystem::parse(Parseur& P) {
	std::string wf(P.get<std::string>("wf"));
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
		ref_(1) = 2;
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
}

void CreateSystem::create(){
	switch(ref_(0)){
		case 2:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:{RGL_ = new ChainFermi(N_,n_,m_,bc_);}break;
									   //case 1:{return ChainDimerized(N_,n_,m_);}break;
								case 2:{RGL_ = new ChainPolymerized(N_,n_,m_,bc_);}break;
								default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						}break;
					default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
				}
			}break;
		case 3:
			//{
			//switch(ref_(1)){
			//case 0:{return TriangleJastrow(N_,n_,m_);}break;
			//case 1:
			//   {
			//   switch(ref_(2)){
			//   case 0:{return TriangleFermi(N_,n_,m_);}break;
			//   case 1:{return TriangleMu(N_,n_,m_);}break;
			//   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
			//   }
			//   }break;
			//case 2:
			//   {
			//   switch(ref_(2)){
			//   case 4:{return TrianglePhi(N_,n_,m_);}break;
			//   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
			//   }
			//   }break;
			//default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
			//}
			//}break;
		case 4:
			{
				switch(ref_(1)){
					//case 0:{return SquareJastrow(N_,n_,m_);}break;
					case 1:
						{
							switch(ref_(2)){
								//   case 0:{return SquareFermi(N_,n_,m_);}break;
								//   case 1:{return SquareMu(N_,n_,m_);}break;
								default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 2:{CGL_ = new SquarePiFlux(N_,n_,m_,bc_);}break;
									   //   case 3:{return SquareSU2PhiFlux(N_,n_,m_);}break;
								default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						}break;
					default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
				}
			}break;
		case 6:
			//{
			//switch(ref_(1)){
			//case 1:
			//{
			//switch(ref_(2)){
			//case 0:{return HoneycombSU3(N_,n_,m_);}break;
			//case 1:{return HoneycombSU4(N_,n_,m_);}break;
			//default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
			//}
			//}break;
			//default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
			//}
			//}break;
		default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
	}
}

void CreateSystem::save(IOFiles& w) const{
	w("ref (type of wavefunction)",ref_);
	if(RGL_){RGL_->save(w);}
	if(CGL_){CGL_->save(w);}
}

void CreateSystem::check(){
	if(RGL_){return RGL_->check();}
	if(CGL_){return CGL_->check();}
}
