#include "CreateSystem.hpp"

CreateSystem::CreateSystem(Parseur& P):
	ref_(3,0),
	N_(P.get<unsigned int>("N")),
	n_(P.get<unsigned int>("n")),
	m_(P.get<unsigned int>("m")),
	M_(N_,(m_*n_ )/N_),
	bc_(P.get<int>("bc")),
	type_(P.get<unsigned int>("type")),
	over_(true),
	RGL_(NULL),
	CGL_(NULL)
{
	parse(P);
}

CreateSystem::~CreateSystem(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
}

void CreateSystem::parse(Parseur& P) {
	std::string wf(P.get<std::string>("wf"));
	if( wf == "chainfermi" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "chainpolymerized" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 1;
		if(P.is_vector("delta")){ 
			Vector<double> tmp(P.get<Vector<double> >("delta"));
			for(unsigned int i(0);i<tmp.size();i++){
				d_.append(tmp(i));
			}
		}
		else { d_.append(P.get<double>("delta")); }
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

	if( wf == "kagomefermi" ){
		ref_(0) = 5;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "kagomedirac" ){
		ref_(0) = 5;
		ref_(1) = 1;
		ref_(2) = 1;
	}
	if( wf == "kagomevbc" ){
		ref_(0) = 5;
		ref_(1) = 2;
		ref_(2) = 0;
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

void CreateSystem::init(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
	std::cout<<"will lunch new system "<<M_<<std::endl;
	switch(ref_(0)){
		case 2:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{RGL_ = new ChainFermi(ref_,N_,m_,n_,M_,bc_);}break;
								case 1:
									{
										RGL_ = new ChainPolymerized(ref_,N_,m_,n_,M_,bc_,d_.last());
										d_.pop();
										if(!d_.size()){ over_ = true; }
									}break;
								default:
									{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						}break;
					default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
				}
			}break;
		case 3:
			//{
			//switch(ref_(1)){
			//case 0:{return TriangleJastrow(N,n,m);}break;
			//case 1:
			//   {
			//   switch(ref_(2)){
			//   case 0:{return TriangleFermi(N,n,m);}break;
			//   case 1:{return TriangleMu(N,n,m);}break;
			//   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
			//   }
			//   }break;
			//case 2:
			//   {
			//   switch(ref_(2)){
			//   case 4:{return TrianglePhi(N,n,m);}break;
			//   default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl; }break;
			//   }
			//   }break;
			//default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
			//}
			//}break;
		case 4:
			{
				switch(ref_(1)){
					//case 0:{return SquareJastrow(N,n,m);}break;
					case 1:
						{
							switch(ref_(2)){
								//   case 0:{return SquareFermi(N,n,m);}break;
								//   case 1:{return SquareMu(N,n,m);}break;
								default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 2:{CGL_ = new SquarePiFlux(ref_,N_,m_,n_,M_,bc_);}break;
								default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						}break;
					default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
				}
			}break;
		case 5:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:{RGL_ = new KagomeFermi(ref_,N_,m_,n_,M_,bc_);}break;
								case 1:{RGL_ = new KagomeDirac(ref_,N_,m_,n_,M_,bc_);}break;
								default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
							}
						} break;
					case 2:
						{
							switch(ref_(2)){
								case 0:{CGL_ = new KagomeVBC(ref_,N_,m_,n_,M_,bc_);}break;
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
			//case 0:{return HoneycombSU3(N,n,m);}break;
			//case 1:{return HoneycombSU4(N,n,m);}break;
			//default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
			//}
			//}break;
			//default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
			//}
			//}break;
		default:{std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;}break;
	}
	switch(type_){
		case 3:
			{
				if(M_(0)){ 
					M_(0) -= 1;
					M_(1) += 1;
					over_ = false;
				} else {
					over_ = true;
				}
			}break;
	}
}
