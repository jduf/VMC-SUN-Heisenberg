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
										RGL_ = new ChainFree(*s_,t,mu);
									}break;
								case 2:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C)    { t = C->get<std::vector<double> >("t"); }
										RGL_ = new ChainPolymerized(*s_,t);
									}break;
								case 3:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C)    { t = C->get<std::vector<double> >("t"); }
										RGL_ = new ChainSAS(*s_,t);
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
									{ RGL_ = new LadderFermi(*s_); }break;
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
								case 2:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderDimerA(*s_,t);
									}break;
								case 3:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderDimerB(*s_,t);
									}break;
								case 4:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderSquarePlaquetteA(*s_,t);
									}break;
								case 5:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderSquarePlaquetteB(*s_,t);
									}break;
								case 6:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderSquarePlaquetteC(*s_,t);
									}break;
								case 7:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderRectangularPlaquetteA(*s_,t);
									}break;
								case 8:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderRectangularPlaquetteB(*s_,t);
									}break;
								case 9:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderRectangularPlaquetteC(*s_,t);
									}break;
								case 10:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new LadderRectangularPlaquetteD(*s_,t);
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
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
								case 4:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new TriangleT3x2(*s_,t);
									}break;
								case 5:
									{
										double t;
										if(param){ t = (*param)(0); }
										if(C){ t = C->get<double>("t"); }
										RGL_ = new TriangleAlternatingPlaquette(*s_,t);
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 1:
									{
										double phi;
										if(param){ phi = (*param)(0); }
										if(C){ phi = C->get<double>("phi"); }
										CGL_ = new TriangleChiral(*s_,phi);
									}break;
								case 2:
									{
										double phi;
										if(param){ phi = (*param)(0); }
										if(C)    { phi = C->get<double>("phi"); }
										CGL_ = new TrianglePhi(*s_,phi);
									}break;
								case 3:
									{
										double phi;
										if(param){ phi = (*param)(0); }
										if(C){ phi = C->get<double>("phi"); }
										CGL_ = new TriangleChiralSG(*s_,phi);
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
							RGL_ = new SquareJastrow(*s_,nu);
						}break;
					case 1:
						{
							switch(ref_(2)){
								case 0:
									{ RGL_ = new SquareFermi(*s_); }break;
								case 1:
									{
										Vector<double> t;
										Vector<double> mu;
										if(param){
											t = (*param);
											mu.set(4,0);
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
								case 3:
									{
										Vector<double> t;
										if(param){ t = (*param); }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareDimerizedBar(*s_,t);
									}break;
								case 4:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareT2x2(*s_,t);
									}break;
								case 5:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareT3x2(*s_,t);
									}break;
								case 6:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareT3x3(*s_,t);
									}break;
								case 7:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareT4x2(*s_,t);
									}break;
								case 8:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareT4x3(*s_,t);
									}break;
								case 9:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareT4x4(*s_,t);
									}break;
								case 10:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareLadder(*s_,t);
									}break;
								case 11:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new SquareVCS(*s_,t);
									}break;
								case 12:
									{
										double mu;
										Vector<double> t;
										if(param){
											mu = (*param)(0);
											t = param->range(1,5);
										}
										if(C){
											mu = C->get<double>("mu");
											t = C->get<std::vector<double> >("t");
										}
										RGL_ = new Squarek2Mu(*s_,mu,t);
									}break;
								case 13:
									{
										double mu;
										Vector<double> t;
										if(param){
											mu = (*param)(0);
											t = param->range(1,5);
										}
										if(C){
											mu = C->get<double>("mu");
											t = C->get<std::vector<double> >("t");
										}
										RGL_ = new Squarek2Mu2x2(*s_,mu,t);
									}break;
								case 14:
									{
										double mu;
										Vector<double> t;
										if(param){
											mu = (*param)(0);
											t = param->range(1,5);
										}
										if(C){
											mu = C->get<double>("mu");
											t = C->get<std::vector<double> >("t");
										}
										RGL_ = new Squarek2Mu2x2Col(*s_,mu,t);
									}break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 1:
									{
										Vector<double> t;
										Vector<double> phi;
										if(param){
											t   = param->range(0,param->size()/2);
											phi = param->range(param->size()/2,param->size());
										}
										if(C){
											t   = C->get<std::vector<double> >("t");
											phi = C->get<std::vector<double> >("phi");
										}
										CGL_ = new SquareFreeFlux(*s_,t,phi);
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
								case 4:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										CGL_ = new SquareBox6(*s_,t);
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
								case 0: { RGL_ = new KagomeFermi(*s_); }break;
								case 1:
										{
											Vector<double> t;
											if(param){ t = (*param); }
											if(C){ t = C->get<std::vector<double> >("t"); }
											RGL_ = new KagomeFree(*s_,t);
										}break;
								case 2:
										{
											double td;
											if(param){ td = (*param)(0); }
											if(C){ td = C->get<double>("td"); }
											RGL_ = new KagomePlaquette3A(*s_,td);
										}break;
								case 3:
										{
											double td;
											if(param){ td = (*param)(0); }
											if(C){ td = C->get<double>("td"); }
											RGL_ = new KagomePlaquette3B(*s_,td);
										}break;
								case 4:
										{
											double td;
											if(param){ td = (*param)(0); }
											if(C){ td = C->get<double>("td"); }
											RGL_ = new KagomePlaquette6A(*s_,td);
										}break;
								case 5:
										{
											double td;
											if(param){ td = (*param)(0); }
											if(C){ td = C->get<double>("td"); }
											RGL_ = new KagomePlaquette6B(*s_,td);
										}break;
								default:{ error(); }break;
							}
						} break;
					case 2:
						{
							switch(ref_(2)){
								case 1:
									{
										double phi;
										if(param){ phi = (*param)(0); }
										if(C){ phi = C->get<double>("phi"); }
										CGL_ = new KagomeChiral(*s_,phi);
									}break;
								case 2: { CGL_ = new KagomeVBC(*s_); }break;
								case 3:
								case 4:
										{
											double phi;
											if(param){ phi = (*param)(0); }
											if(C){ phi = C->get<double>("phi"); }
											CGL_ = new KagomeChiralB(*s_,phi);
										}break;
										{
											Vector<double> t;
											Vector<double> phi;
											if(param){
											}
											if(C){
												t = C->get<std::vector<double> >("t");
												phi = C->get<std::vector<double> >("phi");
											}
											CGL_ = new KagomePiHalfTriangle(*s_,t,phi);
										}break;
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
									{ RGL_ = new HoneycombFermi(*s_); }break;
								case 1:
									{
										Vector<double> t;
										if(param){ t = *param; }
										if(C){ t = C->get<std::vector<double> >("t"); }
										RGL_ = new HoneycombFree(*s_,t);
									}break;
								case 3:
									{
										double t(0);
										unsigned int fc(0);
										if(param){
											t = (*param)(0);
											if(param->size()==2){ fc = (*param)(1); }
											else { fc = 0; }
										}
										if(C) {
											t = C->get<double>("td");
											fc= (C->find("fc",fc,false)?C->get<unsigned int>(fc):0);
										}
										RGL_ = new Honeycomb0pp(*s_,t,fc);
									}break;
								case 4:
									{
										double t(0);
										if(param){ t = (*param)(0); }
										if(C) { t = C->get<double>("td"); }
										RGL_ = new Honeycombp00(*s_,t);
									}break;
								case 5:
									{ RGL_ = new HoneycombPiFlux(*s_); }break;
								default:{ error(); }break;
							}
						}break;
					case 2:
						{
							switch(ref_(2)){
								case 1:
									{
										double phi;
										if(param){ phi = (*param)(0); }
										if(C){ phi = C->get<double>("phi"); }
										CGL_ = new HoneycombChiral(*s_,phi);
									}break;
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
				if(CGL_){
					std::cerr<<__PRETTY_FUNCTION__<<" : might work if C_ knows the parameters"<<std::endl;
					CGL_->create();
				}
			} else { std::cerr<<__PRETTY_FUNCTION__<<" : giving up"<<std::endl; }
		}
	} else {
		if(CGL_){ CGL_->create(); }
		if(CGL_ && get_status()!=1){
			std::cerr<<__PRETTY_FUNCTION__<<" : don't know what to do for the degeneracy at the Fermi level"<<std::endl;
		}
	}
}
/*}*/
