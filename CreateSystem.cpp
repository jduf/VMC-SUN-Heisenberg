#include "CreateSystem.hpp"

CreateSystem::CreateSystem(Parseur& P):
	ref_(3,0),
	N_(P.get<unsigned int>("N")),
	m_(P.get<unsigned int>("m")),
	n_(P.get<unsigned int>("n")),
	M_(N_,(m_*n_ )/N_),
	bc_(P.get<int>("bc")),
	type_(P.get<unsigned int>("type")),
	over_(false),
	RGL_(NULL),
	CGL_(NULL)
{
	parse(P);
}

CreateSystem::CreateSystem(Container const& C):
	ref_(3,0),
	N_(C.get<unsigned int>("N")),
	m_(C.get<unsigned int>("m")),
	n_(C.get<unsigned int>("n")),
	M_(N_,(m_*n_ )/N_),
	bc_(C.get<int>("bc")),
	type_(C.get<unsigned int>("type")),
	over_(false),
	RGL_(NULL),
	CGL_(NULL)
{
	std::string wf(C.get<std::string>("wf"));
	if( wf == "chainpolymerized" ){
		ref_(0) = 2;
		ref_(1) = 1;
		ref_(2) = 1;
		//Vector<double> ti(N_/m_,1);
		//ti(N_/m_-1) = 1-C.get<double>("delta");
		//vd_.add_end(&ti); 
		vd_.add_end(new Vector<double>(N_/m_,1));
		vd_.last()(N_/m_-1) = 1-C.get<double>("delta");
		std::cout<<vd_<<std::endl;
	} else { std::cerr<<"CreateSystem :: lajakdj"<<std::endl; }
}

CreateSystem::CreateSystem(IOFiles* r):
	ref_(r->read<Vector<unsigned int> >()),
	N_(r->read<unsigned int>()),
	m_(r->read<unsigned int>()),
	n_(r->read<unsigned int>()),
	M_(r->read<Vector<unsigned int> >()),
	bc_(r->read<int>()),
	type_(0),
	over_(false),
	RGL_(NULL),
	CGL_(NULL)
{}

CreateSystem::~CreateSystem(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
}

void CreateSystem::parse(Parseur& P){
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
		Vector<double> ti(N_/m_,1);
		if(N_/m_ == 4){
			unsigned int c(0);
			if(P.is_vector("t2")){ c+=1; }
			if(P.is_vector("t4")){ c+=2; }
			switch(c){
				case 0:{
						   ti(1) = P.get<double>("t2");
						   ti(3) = P.get<double>("t4");
						   vd_.add_end(&ti);
					   }break;
				case 1:{
						   Vector<double> tmp(P.get<Vector<double> >("t2"));
						   ti(3) = P.get<double>("t4");
						   for(unsigned int i(0);i<tmp.size();i++){
							   ti(1) = tmp(i);
							   vd_.add_end(&ti);
						   }
					   }break;
				case 2:{
						   Vector<double> tmp(P.get<Vector<double> >("t4"));
						   ti(1) = P.get<double>("t2");
						   for(unsigned int i(0);i<tmp.size();i++){
							   ti(3) = tmp(i);
							   vd_.add_end(&ti);
						   }
					   }break;
				default:{std::cerr<<"void CreateSystem::parse(Parseur& P) : unknown tij"<<std::endl; }
			}
		} else { 
			if(P.is_vector("delta")){ 
				Vector<double> tmp(P.get<Vector<double> >("delta"));
				for(unsigned int i(0);i<tmp.size();i++){
					ti(N_/m_-1) = 1-tmp(i);
					vd_.add_end(&ti);
				}
			} else {
				ti(N_/m_-1) = 1-P.get<double>("delta");
				vd_.add_end(&ti); 
			}
		}
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
		sel0_ = 0;
		sel1_ = 0;
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
		if(P.is_vector("td")){ 
			Vector<double> tmp(P.get<Vector<double> >("td"));
			for(unsigned int i(0);i<tmp.size();i++){
				d_.add_end(&tmp(i));
			}
		}
		else { d_.add_end(new double(P.get<double>("td"))); }
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
		case 2:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{
										RGL_ = new ChainFermi<double>(ref_,N_,m_,n_,M_,bc_);
										over_ = true;
									}break;
								case 1:
									{
										if(read){ RGL_ = new ChainPolymerized(ref_,N_,m_,n_,M_,bc_,read->read<Vector<double> >()); } 
										else { 
											RGL_ = new ChainPolymerized(ref_,N_,m_,n_,M_,bc_,vd_.last()); 
											vd_.pop_end();
											if(!vd_.size()){ over_ = true; }
										}
									}break;
								default: {error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{
										CGL_ = new ChainFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_);
										over_ = true;
									}break;
								default: {error();}break;
							}
						}break;
					default:{error();}break;
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
									   //   case 1:{return TriangleMu(N,n,m);}break;
								default:{error();}break;
							}
						}break;
						//case 2:
						//   {
						//   switch(ref_(2)){
						//   case 4:{return TrianglePhi(N,n,m);}break;
						//   default:{error();}break;
						//   }
						//   }break;
						//default:{error();}break;
				}
			}break;
		case 4:
			{
				switch(ref_(1)){
					case 0:
						{RGL_ = new SquareJastrow(ref_,N_,m_,n_,M_,bc_);}break;
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{
										RGL_ = new SquareFermi<double>(ref_,N_,m_,n_,M_,bc_);
										over_ = true;
									}break;
									//   case 1:{return SquareMu(N,n,m);}break;
								default:{error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{
										CGL_ = new SquareFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_);
										over_ = true;
									}break;
								case 2:{CGL_ = new SquarePiFlux(ref_,N_,m_,n_,M_,bc_);}break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		case 5:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:{RGL_ = new KagomeFermi<double>(ref_,N_,m_,n_,M_,bc_,sel0_,sel1_);}break;
								case 1:{RGL_ = new KagomeDirac<double>(ref_,N_,m_,n_,M_,bc_);}break;
								default:{error();}break;
							}
						} break;
					case 2:
						{
							switch(ref_(2)){
								case 0:{CGL_ = new KagomeFermi<std::complex<double> >(ref_,N_,m_,n_,M_,bc_,sel0_,sel1_);}break;
								case 1:{CGL_ = new KagomeDirac<std::complex<double> >(ref_,N_,m_,n_,M_,bc_);}break;
								case 2:{CGL_ = new KagomeVBC(ref_,N_,m_,n_,M_,bc_);}break;
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
									{
										if(read){ d_.add_end( new double(read->read<double>()) ); }
										RGL_ = new Honeycomb0pp(ref_,N_,m_,n_,M_,bc_,d_.last());
										d_.pop_end();
										if(!d_.size()){ over_ = true; }
									}break;
									//case 1:{return HoneycombSU4(N,n,m);}break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		default:{error();}break;
	}
	switch(type_){
		case 2:
			{
				if(CGL_){
					sel0_++;
					over_ = false;
					if(sel0_==9){ 
						sel1_++;
						sel0_ = 0;
						if(sel1_==639){ over_ = true; }
						//sel0_ = sel1_;
						//if(sel1_==9){ over_ = true; }
					}
				}
			}break;
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
	if(ios){
		if(RGL_){ RGL_->set_IOSystem(ios); }
		if(CGL_){ CGL_->set_IOSystem(ios); }
	}
}

void CreateSystem::error(){
	std::cerr<<"ref_ = ["<<ref_(0)<<ref_(1)<<ref_(2)<<"] unknown"<<std::endl;
	over_ = true;
}  

void CreateSystem::create(){
	if(RGL_){
		RGL_->create();
		if(get_status()!=1 && type_ == 1){
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
