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
	if(RGL_){ delete RGL_; }
	if(CGL_){ delete CGL_; }
}
/*}*/

/*core methods*/
/*{*/
void CreateSystem::init(Vector<double> const* const param, Container* C){
	if(RGL_){ delete RGL_; }
	if(CGL_){ delete CGL_; }
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
										Vector<double> t;
										Vector<double> mu;
										if(param){
											t.set(param->size()/2);
											mu.set(param->size()/2);
											for(unsigned int i(0);i<t.size();i++){ t(i) = (*param)(i); }
											for(unsigned int i(0);i<mu.size();i++){ mu(i) = (*param)(i+t.size()); }
										}
										if(C){
											t = C->get<std::vector<double> >("t");
											mu = C->get<std::vector<double> >("mu");
										}
										if(t.ptr()){ RGL_ = new ChainFree(*s_,t,mu); }
									}break;
								case 2:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C)    { t = C->get<std::vector<double> >("t"); }
										if(t.ptr()){ RGL_ = new ChainPolymerized(*s_,t); }
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0: 
									{ CGL_ = new ChainFermi<std::complex<double> >(*s_); }break;
								default:{ error(); }break;
							}
						}break;
					default:{ error(); }break;
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
								case 1:
									{
										Vector<double> t;
										Vector<double> mu;
										if(param){
											mu.set((param->size()+1)*2/5);
											t.set(param->size()-mu.size());
											for(unsigned int i(0);i<t.size();i++){ t(i) = (*param)(i); }
											for(unsigned int i(0);i<mu.size();i++){ mu(i) = (*param)(i+t.size()); }
										}
										if(C){
											t = C->get<std::vector<double> >("t");
											mu = C->get<std::vector<double> >("mu");
										}
										RGL_ = new LadderFree(*s_,t,mu);
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{ CGL_ = new LadderFermi<std::complex<double> >(*s_); }break;
								case 1:
									{
										Vector<double> t;
										Vector<double> flux;
										if(param){
											unsigned int l((3*param->size()-1)/4);
											t.set(l);
											flux.set(param->size()-l);
											for(unsigned int i(0);i<l;i++){ t(i) = (*param)(i); }
											for(unsigned int i(0);i<param->size()-l;i++){ flux(i) = (*param)(l+i); }
										}
										if(C){
											t = C->get<std::vector<double> >("t");
											flux = C->get<std::vector<double> >("flux");
										}
										CGL_ = new LadderFreeFlux(*s_,t,flux);
									}break;
								default:{ error(); }break;
							}
						}break;
					default:{ error(); }break;
				}
			}break;
		case 3:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ RGL_ = new TriangleFermi(*s_); }break;
								case 1:
									{
										Vector<double> t;
										Vector<double> mu;
										if(param){
										}
										if(C){
											t = C->get<std::vector<double> >("t");
											mu = C->get<std::vector<double> >("mu");
										}
										RGL_ = new TriangleFree(*s_,t,mu);
									}break;
								case 2:
									{
										double t;
										if(param){ t = (*param)(0); }
										if(C){ t = C->get<double>("t"); }
										RGL_ = new TrianglePlaquette(*s_,t);
									}break;
								case 3:
									{
										double mu;
										if(param){ mu = (*param)(0); }
										if(C){ mu = C->get<double>("mu"); }
										RGL_ = new TriangleMu(*s_,mu);
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 1:
									{ CGL_ = new TriangleChiral(*s_); }break;
								case 2:
									{
										double phi;
										if(param){ phi = (*param)(0); }
										if(C)    { phi = C->get<double>("phi"); }
										CGL_ = new TrianglePhi(*s_,phi);
									}break;
								default:{ error(); }break;
							}
						}break;
					default:{ error(); }break;
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
								case 1:
									{
										Vector<double> t(2);
										Vector<double> mu(1);
										if(param){
											mu(0)= (*param)(0);
											t(0) = (*param)(1);
											t(1) = (*param)(2);
										}
										if(C){
											t = C->get<std::vector<double> >("t");
											mu = C->get<std::vector<double> >("mu");
										}
										RGL_ = new SquareFree(*s_,t,mu);
									}break;
								case 2:
									{
										double mu;
										if(param){ mu = (*param)(0); }
										if(C){ mu = C->get<double>("mu"); }
										RGL_ = new SquareMu(*s_,mu);
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 1:
									{
										Vector<double> phi;
										if(param){ phi = param->range(0,0); }
										if(C){ phi=C->get<std::vector<double> >("phi"); }
										CGL_ = new SquareFreeFlux(*s_,phi);
									}break;
								case 2:
									{ CGL_ = new SquarePiFlux(*s_); }break;
								case 3:
									{ 
										double phi;
										if(param){ phi = (*param)(0); }
										if(C){ phi = C->get<double>("phi"); }
										CGL_ = new SquareChiral(*s_,phi);
									}break;
								default:{ error(); }break;
							}
						}break;
					default:{ error(); }break;
				}
			}break;
		case 5:
			{
				switch(ref_(1)){
					case 1:
						{
							switch(ref_(2)){
								case 0: { RGL_ = new KagomeFermi<double>(*s_); }break;
								case 1: { RGL_ = new KagomeDirac<double>(*s_); }break;
								default:{ error(); }break;
							}
						} break;
					case 2:
						{
							switch(ref_(2)){
								case 0: { CGL_ = new KagomeFermi<std::complex<double> >(*s_); }break;
								case 1: { CGL_ = new KagomeDirac<std::complex<double> >(*s_); }break;
								case 2: { CGL_ = new KagomeVBC(*s_); }break;
								default:{ error(); }break;
							}
						}break;
					default:{ error(); }break;
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
										double t(0);
										if(param){ t = (*param)(0); }
										if(C)    { t = C->get<double>("td"); }
										if(param || C){ RGL_ = new Honeycomb0pp(*s_,t); }
									}break;
								case 1:
									{ RGL_ = new HoneycombPiFlux(*s_); }break;
								case 2:
									{ 
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new HoneycombFree(*s_,t); 
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 0:
									{ CGL_ = new HoneycombChiral(*s_); }break;
								default:{ error(); }break;
							}
						}break;
					default:{ error(); }break;
				}
			}break;
		default:{ error(); }break;
	}
	//!Handle geometry problem
	if(RGL_ && RGL_->get_status()==3){
		delete RGL_;
		RGL_ = NULL;
		ref_(3)++;
		if(s_->try_other_geometry(ref_)){ init(param,C); }
	}
	if(CGL_ && CGL_->get_status()==3){
		delete CGL_;
		CGL_ = NULL;
		ref_(3)++;
		if(s_->try_other_geometry(ref_)){ init(param,C); }
	}
}

void CreateSystem::create(bool const& try_solve_degeneracy){
	if(RGL_){
		RGL_->create();
		if(RGL_->get_status()!=1){
			if(try_solve_degeneracy){
				delete RGL_;
				RGL_ = NULL;
				ref_(1)=2;
				init(NULL,C_);
				std::cerr<<__PRETTY_FUNCTION__<<" : would work if C_ knows the parameters"<<std::endl;
				if(CGL_){ CGL_->create(); }
			} else { std::cerr<<__PRETTY_FUNCTION__<<" : giving up"<<std::endl; }
		}
	} else {
		if(CGL_){ CGL_->create(); }
		if(CGL_ && get_status()!=1){
			std::cerr<<__PRETTY_FUNCTION__<<" : behaviour undefined"<<std::endl;
		}
	}
}
/*}*/
