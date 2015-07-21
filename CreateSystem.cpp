#include "CreateSystem.hpp"

CreateSystem::CreateSystem(System* s):
	s_(s),
	ref_(s_->get_ref()),
	RGL_(NULL),
	CGL_(NULL)
{}

CreateSystem::~CreateSystem(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
}

Vector<unsigned int> CreateSystem::get_ref(std::string const& wf){
	Vector<unsigned int> ref(3,6);
	if( wf == "chainfermi" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "chainpolymerized" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 1;
	}

	if( wf == "ladderfermi"){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "ladderfree"){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 4;
	}

	if( wf == "trianglefermi" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "trianglemu" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "trianglephi" ){
		ref(0) = 3;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "trianglejastrow" ){
		ref(0) = 3;
		ref(1) = 0;
	}

	if( wf == "squarefermi" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "squaremu" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 1;
	}
	//if( wf == "squarefreereal" ){
	//ref(0) = 4;
	//ref(1) = 1;
	//ref(2) = 3;
	//C_.set("t",param?param->range(0,3):C->get<std::vector<double> >("t"));
	//C_.set("mu",param?param->range(3,5):C->get<std::vector<double> >("mu"));
	//}
	if( wf == "squareacsl" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
	}
	if( wf == "squarefreecomplex" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 4;
	}
	if( wf == "squarepiflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "squarephi" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
	}
	if( wf == "squarejastrow" ){
		ref(0) = 4;
		ref(1) = 0;
	}

	if( wf == "kagomefermi" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "kagomedirac" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "kagomevbc" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 2;
	}

	if( wf == "honeycomb0pp" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "honeycombsu4" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 0;
	}
	return ref;
}

void CreateSystem::set_param(Container* C, Vector<double> const* param){
	switch(ref_(0)){
		case 1:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 1:
									{ C_.set("t",param?(*param):C->get<double>("t")); }break;
								default: {error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		case 2:	
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 4:
									{ C_.set("t",param?(*param):C->get<std::vector<double> >("t")); }break;
								default: {error();}break;
							}
						}break;
					default: {error();}break;
				}
			}break;
		case 4:
			{
				switch(ref_(1)){
					case 0:
						{
							Vector<double> tmp(C->get<Vector<double> >("nu"));
							Matrix<double> nu(4,2);
							for(unsigned int i(0);i<4;i++){
								nu(i,0) = tmp(0);
								nu(i,1) = tmp(1);
							}
							C_.set("nu",nu);
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 3:
									{ C_.set("t",param?*param:C->get<std::vector<double> >("t")); }break;
								case 4:
									{
										/*3 free parameters*/
										C_.set("mu",param?param->range(0,0):C->get<std::vector<double> >("mu"));
										C_.set("t",param?param->range(1,2):C->get<std::vector<double> >("t"));
										C_.set("phi",Vector<double>(1,M_PI/4.0));

										/*1 free parameter, box*/
										//C_.set("mu",Vector<double>(1,0));
										//C_.set("t",param?Vector<double>(2,(*param)(0)):C->get<std::vector<double> >("t"));
										//C_.set("phi",Vector<double>(1,M_PI/4.0));

										/*1 free parameter, dimer*/
										//Vector<double> t(2);
										//t(0) = (*param)(0);
										//t(1) = 1;
										//C_.set("mu",Vector<double>(1,0));
										//C_.set("t",param?t:C->get<std::vector<double> >("t"));
										//C_.set("phi",Vector<double>(1,M_PI/4.0));
									}break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		case 6:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ C_.set("td",param?(*param)(0):C->get<double>("td")); }break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		default:{error();}break;
	}
}

void CreateSystem::init(IOFiles* read, IOSystem* ios){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
	switch(ref_(0)){
		case 1:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ RGL_ = new ChainFermi<double>(*s_); }break;
								case 1:
									{
										if(read){ RGL_ = new ChainPolymerized(*s_,read->read<Vector<double> >()); }
										else    { RGL_ = new ChainPolymerized(*s_,C_.get<Vector<double> >("t")); }
									}break;
								default: {error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{ CGL_ = new ChainFermi<std::complex<double> >(*s_); }break;
								default: {error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		case 2:	
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ RGL_ = new LadderFermi<double>(*s_); }break;
								case 4:
									{ RGL_ = new LadderFree(*s_,C_.get<Vector<double> >("t")); }break;
								default: {error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{ CGL_ = new LadderFermi<std::complex<double> >(*s_); }break;
								default: {error();}break;
							}
						}break;
					default: {error();}break;
				}
			}break;
		case 3:
			{
				switch(ref_(1)){
					//case 0:{return TriangleJastrow(N,n,m);}break;
					case 1:
						{
							switch(ref_(2)){
								case 0:{RGL_ = new TriangleFermi(*s_);}break;
									   //   //   case 1:{return TriangleMu(N,n,m);}break;
								default:{error();}break;
							}
						}break;
						////case 2:
						////   {
						////   switch(ref_(2)){
						////   case 4:{return TrianglePhi(N,n,m);}break;
						////   default:{error();}break;
						////   }
						////   }break;
						////default:{error();}break;
				}
			}break;
		case 4:
			{
				switch(ref_(1)){
					case 0:
						{RGL_ = new SquareJastrow(*s_,C_.get<Matrix<double> >("nu"));}break;
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ RGL_ = new SquareFermi<double>(*s_); }break;
								default:{error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 2:
									{ CGL_ = new SquarePiFlux(*s_); }break;
								case 3:
									{ CGL_ = new SquareACSL(*s_,C_.get<Vector<double> >("t")); }break;
								case 4:
									{ CGL_ = new SquareFreeComplex(*s_,C_.get<Vector<double> >("t"),C_.get<Vector<double> >("mu"),C_.get<Vector<double> >("phi")); }break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
			//case 5:
			//{
			//switch(ref_(1)){
			//case 1:
			//{
			//switch(ref_(2)){
			//case 0:{
			//   std::cerr<<"KagomeFermi<double>(*s_,Vector<unsigned int>,Vector<unsigned int>) not fully defined"<<std::endl;
			//   //   RGL_ = new KagomeFermi<double>(*s_,sel0_,sel1_);
			//   }break;
			//case 1:{RGL_ = new KagomeDirac<double>(*s_);}break;
			//default:{error();}break;
			//}
			//} break;
			//case 2:
			//{
			//switch(ref_(2)){
			//case 0:{
			//   std::cerr<<"KagomeFermi<std::complex<double> >(*s_,Vector<unsigned int>,Vector<unsigned int>) not fully defined"<<std::endl;
			//   //   CGL_ = new KagomeFermi<std::complex<double> >(*s_,sel0_,sel1_);
			//   }break;
			//case 1:{CGL_ = new KagomeDirac<std::complex<double> >(*s_);}break;
			//case 2:{CGL_ = new KagomeVBC(ref_,N_,m_,n_,M_,bc_);}break;
			//default:{error();}break;
			//}
			//}break;
			//default:{error();}break;
			//}
			//}break;
		case 6:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{
										if(read){ RGL_ = new Honeycomb0pp(*s_,read->read<double>()) ; }
										else    { RGL_ = new Honeycomb0pp(*s_,C_.get<double>("td")); }
									}break;
									////case 1:{return HoneycombSU4(N,n,m);}break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		default:{error();}break;
	}
	if(ios){
		if(RGL_){ RGL_->set_IOSystem(ios); }
		if(CGL_){ CGL_->set_IOSystem(ios); }
	}
}

void CreateSystem::error() const {
	std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;
}

void CreateSystem::create(bool try_solve_degeneracy){
	if(RGL_){
		RGL_->create();
		if(get_status()!=1){
			if(try_solve_degeneracy){
				delete RGL_;
				RGL_=NULL;
				ref_(1)=2;
				init();
				if(CGL_){ CGL_->create(); }
			} else {
				std::cerr<<"void CreateSystem::create(bool try_solve_degeneracy) : giving up"<<std::endl;
			}
		}
	} else {
		if(CGL_){ CGL_->create(); }
		if(CGL_ && get_status()!=1){
			std::cerr<<"void CreateSystem::create(bool try_solve_degeneracy) : behaviour undefined"<<std::endl;
		}
	}
}
