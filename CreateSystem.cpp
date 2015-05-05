#include "CreateSystem.hpp"

CreateSystem::CreateSystem(Container* C, Vector<double> const* param):
	ref_(3,0),
	N_(C->get<unsigned int>("N")),
	m_(C->get<unsigned int>("m")),
	n_(C->get<unsigned int>("n")),
	M_(N_,(m_*n_ )/N_),
	bc_(C->get<int>("bc")),
	RGL_(NULL),
	CGL_(NULL)
{
	parse(C,param);
}

CreateSystem::CreateSystem(IOFiles* r):
	ref_(*r),
	N_(r->read<unsigned int>()),
	m_(r->read<unsigned int>()),
	n_(r->read<unsigned int>()),
	M_(*r),
	bc_(r->read<int>()),
	RGL_(NULL),
	CGL_(NULL)
{}

CreateSystem::~CreateSystem(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
}

void CreateSystem::parse(Container* C, Vector<double> const* param){
	std::string wf(C->get<std::string>("wf"));
	if( wf == "chainfermi" ){
		ref_(0) = 1;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "chainpolymerized" ){
		ref_(0) = 1;
		ref_(1) = 1;
		ref_(2) = 1;
		Vector<double> t(N_/m_,1);
		if(N_/m_ == 4){
			t(1) = param?(*param)(0):C->get<double>("t2");
			t(3) = param?(*param)(1):C->get<double>("t4");
		} else { 
			t(N_/m_-1) = param?(*param)(0):C->get<double>("t2");
		}
		C_.set("t",t);
	}
	if( wf == "chainpolymerizedjjp" ){
		if(N_/m_ == 2){
			ref_(0) = 2;
			ref_(1) = 1;
			ref_(2) = 2;
			Vector<double> t(N_/m_,0);
			t(0) = C->get<double>("t");
			t(1) = C->get<double>("tp");
			Vector<double> J(2);
			J(0) = C->get<double>("J");
			J(1) = C->get<double>("Jp");

			C_.set("t",t);
			C_.set("J",J);
		} else { 
			std::cerr<<"chainpolymerizedjjp N/m!=2"<<std::endl;
		}
	}

	if( wf == "ladderfermi"){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 0;
	}
	if( wf == "ladderfree"){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 4;
		C_.set("t",param?(*param):C->get<std::vector<double> >("t"));
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
	//if( wf == "squarefreereal" ){
	//ref_(0) = 4;
	//ref_(1) = 1;
	//ref_(2) = 3;
	//C_.set("t",param?param->range(0,3):C->get<std::vector<double> >("t"));
	//C_.set("mu",param?param->range(3,5):C->get<std::vector<double> >("mu"));
	//}
	if( wf == "squarefreecomplex" ){
		ref_(0) = 4;
		ref_(1) = 2;
		ref_(2) = 3;

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
	}
	if( wf == "squarepiflux" ){
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

		Vector<double> tmp(C->get<Vector<double> >("nu"));
		Matrix<double> nu(4,2);
		for(unsigned int i(0);i<4;i++){
			nu(i,0) = tmp(0);
			nu(i,1) = tmp(1);
		}
		C_.set("nu",nu);
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
		ref_(2) = 2;
	}

	if( wf == "honeycomb0pp" ){
		ref_(0) = 6;
		ref_(1) = 1;
		ref_(2) = 0;
		C_.set("td",param?(*param)(0):C->get<double>("td"));
	}
	if( wf == "honeycombsu4" ){
		ref_(0) = 6;
		ref_(1) = 1;
		ref_(2) = 0;
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
									{ RGL_ = new ChainFermi<double>(ref_,N_,m_,n_,M_,bc_); }break;
								case 1:
									{
										if(read){ RGL_ = new ChainPolymerized(ref_,N_,m_,n_,M_,bc_,read->read<Vector<double> >()); } 
										else { RGL_ = new ChainPolymerized(ref_,N_,m_,n_,M_,bc_,C_.get<Vector<double> >("t")); }
									}break;
								case 2:
									{
										if(read){ RGL_ = new ChainPolymerizedJJp(ref_,N_,m_,n_,M_,bc_,read->read<Vector<double> >(),read->read<Vector<double> >()); } 
										else { RGL_ = new ChainPolymerizedJJp(ref_,N_,m_,n_,M_,bc_,C_.get<Vector<double> >("t"),C_.get<Vector<double> >("J")); }
									}break;
								default: {error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{ CGL_ = new ChainFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_); }break;
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
									{ RGL_ = new LadderFermi<double>(ref_,N_,m_,n_,M_,bc_); }break;
								case 4:
									{ RGL_ = new LadderFree(ref_,N_,m_,n_,M_,bc_,C_.get<Vector<double> >("t")); }break;
								default: {error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{ CGL_ = new LadderFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_); }break;
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
								case 0:{RGL_ = new TriangleFermi(ref_,N_,m_,n_,M_,bc_);}break;
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
						{RGL_ = new SquareJastrow(ref_,N_,m_,n_,M_,bc_,C_.get<Matrix<double> >("nu"));}break;
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ RGL_ = new SquareFermi<double>(ref_,N_,m_,n_,M_,bc_); }break;
									//case 3:
									//{ RGL_ = new SquareFreeReal(ref_,N_,m_,n_,M_,bc_,C_.get<Vector<double> >("t"),C_.get<Vector<double> >("mu")); }break;
									//   case 1:{return SquareMu(N,n,m);}break;
								default:{error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								//case 0:
								//{ CGL_ = new SquareFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_); }break;
								case 2:
									{ CGL_ = new SquarePiFlux(ref_,N_,m_,n_,M_,bc_); }break;
								case 3:
									{ CGL_ = new SquareFreeComplex(ref_,N_,m_,n_,M_,bc_,C_.get<Vector<double> >("t"),C_.get<Vector<double> >("mu"),C_.get<Vector<double> >("phi")); }break;
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
			//   std::cerr<<"KagomeFermi<double>(ref_,N_,m_,n_,M_,bc_,Vector<unsigned int>,Vector<unsigned int>) not fully defined"<<std::endl;
			//   //   RGL_ = new KagomeFermi<double>(ref_,N_,m_,n_,M_,bc_,sel0_,sel1_);
			//   }break;
			//case 1:{RGL_ = new KagomeDirac<double>(ref_,N_,m_,n_,M_,bc_);}break;
			//default:{error();}break;
			//}
			//} break;
			//case 2:
			//{
			//switch(ref_(2)){
			//case 0:{
			//   std::cerr<<"KagomeFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_,Vector<unsigned int>,Vector<unsigned int>) not fully defined"<<std::endl;
			//   //   CGL_ = new KagomeFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_,sel0_,sel1_);
			//   }break;
			//case 1:{CGL_ = new KagomeDirac<std::complex<double> >(ref_,N_,m_,n_,M_,bc_);}break;
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
										if(read){ RGL_ = new Honeycomb0pp(ref_,N_,m_,n_,M_,bc_,read->read<double>()) ; }
										else { RGL_ = new Honeycomb0pp(ref_,N_,m_,n_,M_,bc_,C_.get<double>("td")); }
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

void CreateSystem::create(){
	if(RGL_){
		RGL_->create();
		if(get_status()!=1){
			delete RGL_;
			RGL_=NULL;
			ref_(1)=2;
			init();
			if(CGL_){ CGL_->create(); }
		}
	} else {
		if(CGL_){ CGL_->create(); }
		if(CGL_ && get_status()!=1){
			std::cerr<<"void CreateSystem::create() : behaviour undefined"<<std::endl;
		}
	}
}  
