#include "CreateSystem.hpp"

/*constructors, destructors*/
/*{*/
CreateSystem::CreateSystem(System const* const s):
	s_(s),
	ref_(s_->get_ref()),
	RGL_(NULL),
	CGL_(NULL)
{}

CreateSystem::~CreateSystem(){
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
}
/*}*/

/*core methods*/
/*{*/
void CreateSystem::init(Vector<double> const* const param, Container* C){
	if(C){ C_ = C; }
	if(RGL_){delete RGL_;}
	if(CGL_){delete CGL_;}
	switch(ref_(0)){
		case 1:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0: { RGL_ = new ChainFermi<double>(*s_); }break;
								case 1:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C)    { t = C->get<std::vector<double> >("t"); }
										if(t.ptr()){ RGL_ = new ChainPolymerized(*s_,t); }
									}break;
								default:{error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0: { CGL_ = new ChainFermi<std::complex<double> >(*s_); }break;
								default:{error();}break;
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
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C)    { t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderFree(*s_,t);
									}break;
								default:{error();}break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0: { CGL_ = new LadderFermi<std::complex<double> >(*s_); }break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		case 3:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0: { RGL_ = new TriangleFermi(*s_); }break;
								default:{error();}break;
							}
						}break;
					default:{error();}break;
				}
			}break;
		case 4:
			{
				switch(ref_(1)){
					case 0:
						{
							Vector<double> tmp;
							Matrix<double> nu;
							if(C){
								tmp = C->get<std::vector<double> >("nu");
								nu.set(4,2);
								for(unsigned int i(0);i<4;i++){
									nu(i,0) = tmp(0);
									nu(i,1) = tmp(1);
								}
							}
							if(nu.ptr()){ RGL_ = new SquareJastrow(*s_,nu); }
						}break;
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
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C)    { t = C->get<std::vector<double> >("t"); }
										if(t.ptr()){ CGL_ = new SquareACSL(*s_,t); }
									}break;
								case 4:
									{
										Vector<double> t;
										Vector<double> mu;
										Vector<double> phi;
										if(param){
											t = param->range(0,0);
											mu= param->range(1,2);
											phi.set(1,M_PI/4.0);
										}
										if(C){
											t = C->get<std::vector<double> >("t");
											mu= C->get<std::vector<double> >("mu");
											phi=Vector<double>(1,M_PI/4.0);
										}
										if(mu.ptr() && mu.ptr() && t.ptr()){
											CGL_ = new SquareFreeComplex(*s_,t,mu,phi);
										}
									}break;
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
								case 1: { RGL_ = new KagomeDirac<double>(*s_); }break;
								default:{error();}break;
							}
						} break;
					case 2:
						{
							switch(ref_(2)){
								case 1: {CGL_ = new KagomeDirac<std::complex<double> >(*s_); }break;
								case 2: {CGL_ = new KagomeVBC(*s_); }break;
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
										double t;
										if(param){ t = (*param)(0); }
										if(C)    { t = C->get<double>("td"); }
										if(param || C){ RGL_ = new Honeycomb0pp(*s_,t); }
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
}

void CreateSystem::create(bool const& try_solve_degeneracy){
	if(RGL_){
		RGL_->create();
		if(RGL_->get_status()!=1){
			if(try_solve_degeneracy){
				ref_(1)=2;
				init(NULL,C_);
				std::cerr<<__PRETTY_FUNCTION__<<" : need to check if this works"<<std::endl;
				if(CGL_){ CGL_->create(); }
			} else {
				std::cerr<<__PRETTY_FUNCTION__<<" : giving up"<<std::endl;
			}
		}
	} else {
		if(CGL_){ CGL_->create(); }
		if(CGL_ && get_status()!=1){
			std::cerr<<__PRETTY_FUNCTION__<<" : behaviour undefined"<<std::endl;
		}
	}
}
/*}*/
