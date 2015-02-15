#include"PSO.hpp"

/*constructors and destructor*/
/*{*/
PSO::PSO(unsigned int Nparticles, unsigned int Nfreedom, double cg, double cp, unsigned int maxiter):
	Nparticles_(Nparticles),
	bparticle_(0),
	Nfreedom_(Nfreedom),
	maxiter_(maxiter),
	min_(Nfreedom,-1.5),
	max_(Nfreedom,1.5),
	pbx_(new Vector<double>[Nparticles]),
	pv_(new Vector<double>[Nparticles]),
	px_(new Vector<double>[Nparticles]),
	free_(new bool[Nparticles]),
	pfbx_(new double[Nparticles]),
	cg_(cg),
	cp_(cp),
	chi_(-2.0/(2.0-(cp+cg)-sqrt((cp+cg)*(cp+cg)-4*(cp+cg)))),
	rnd_(0.0,1.0)
{
	if(cg_+cp_<4){
		std::cerr<<"PSO::PSO(Nparam,cg,cp,*f) : cg+cp<4 => chi=nan. Redefinition cg=cp=2.05"<<std::endl;
		cg_ = 2.05;
		cp_ = 2.05;
	}
	for(unsigned int i(0);i<Nparticles_;i++){
		px_[i].set(Nfreedom_);
		pv_[i].set(Nfreedom_);
		pbx_[i].set(Nfreedom_);
		free_[i] = true;
	}
}

void PSO::PSO_init(){
#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int i=0;i<Nparticles_;i++){
#pragma omp critical
		{
			for(unsigned int j(0);j<Nfreedom_;j++){
				px_[i](j) = rnd_.get()*(max_(j)-min_(j))+min_(j);
				pv_[i](j) = rnd_.get()*(max_(j)-min_(j))+min_(j);
			}
		}
		pbx_[i] = px_[i];
		pfbx_[i] = f(px_[i]);
		std::cout<<px_[i]<<std::endl;
	}
	for(unsigned int i(1);i<Nparticles_;i++){
		if(pfbx_[i] < pfbx_[bparticle_] ){
			bparticle_ = i;
		}
	}
}

void PSO::PSO_set_limit(unsigned int param, double min, double max){
	min_(param) = min;
	max_(param) = max;
}

PSO::~PSO(){
	if(pbx_){delete[] pbx_;}
	if(px_){delete[] px_;}
	if(pv_){delete[] pv_;}
	if(pfbx_){delete[] pfbx_;}
	if(free_){delete[] free_;}
}
/*}*/

/*core of the class*/
/*{*/
void PSO::move(unsigned int i){
	for(unsigned int j(0);j<Nfreedom_;j++){
		pv_[i](j) = chi_*(pv_[i](j) + cp_*rnd_.get()*(pbx_[i](j)-px_[i](j)) + cg_*rnd_.get()*(pbx_[bparticle_](j)-px_[i](j)));
		//if(v_[i] > (max_[i]-min_[i])/2.0){v_[i] = (max_[i]-min_[i])/4.0;}
		//if(v_[i] < (min_[i]-max_[i])/2.0){v_[i] = (min_[i]-max_[i])/4.0;}
		if( px_[i](j)+pv_[i](j) > max_(j)){ 
			pv_[i](j) = log(1.0+rnd_.get()*(exp(max_(j)-px_[i](j))-1.0));
		}
		if( px_[i](j)+pv_[i](j) < min_(j)){ 
			pv_[i](j) =-log(1.0+rnd_.get()*(exp(px_[i](j)-min_(j))-1.0));
		}
		px_[i](j) += pv_[i](j); 
	}
}

void PSO::move_on_grid(unsigned int i){
	unsigned int n;
	double dx(0.01);
	for(unsigned int j(0);j<Nfreedom_;j++){
		n=0;
		if(std::abs(px_[i](j))<dx/2){ n=1; }
		if(std::abs(px_[i](j)-min_(j))<dx/2){ n=2; }
		if(std::abs(px_[i](j)-max_(j))<dx/2){ n=3; }
		switch(n){
			case 0:{ px_[i](j) = std::round(std::abs(px_[i](j)/dx))*dx; }break;
			case 1:{ px_[i](j) = 0; }break;
			case 2:{ px_[i](j) = min_(j); }break;
			case 3:{ px_[i](j) = max_(j); }break;
		}
	}
}

void PSO::evaluate(unsigned int i){
	double fx(f(px_[i]));
	if( fx < pfbx_[i]){
		pfbx_[i] = fx; 
		pbx_[i] = px_[i]; 
	}
}

void PSO::next_step(unsigned int i){
#pragma omp critical
	{
		move(i);
	}
	move_on_grid(i);
	evaluate(i);
	if(pfbx_[i] < pfbx_[bparticle_] ){
		bparticle_ = i;
	}
}

void PSO::PSO_run(bool synchro){
	std::cout<<"PSO::will start the run"<<std::endl;	
	if(synchro){
		for(unsigned int iter(0); iter<maxiter_; iter++){
#pragma omp parallel for
			for(unsigned int i=0;i<Nparticles_;i++){
				next_step(i);
			}
		}
	} else {
		//std::cerr<<"Nparticles must be bigger than omp_get_num_threads"<<std::endl;
		unsigned int i(0);
		unsigned int ip(0);
#pragma omp parallel for schedule(dynamic,1) firstprivate(ip)
		for(unsigned int iter=0; iter<maxiter_; iter++){
#pragma omp critical
			{
				i++;
				ip = (i-1) % Nparticles_;
				if(!free_[ip]){ iter--; }
			}
			if(free_[ip]){
				free_[ip]=false;
				next_step(ip);
				free_[ip]=true;
			}
		}
	}
}
/*}*/

/*print, save and load*/
/*{*/
void PSO::PSO_print(){
	for(unsigned int i(0);i<Nparticles_;i++){
		std::cout<<px_[i]<<pv_[i]<<pbx_[i]<<pfbx_[i]<<std::endl;
	}
}

void PSO::PSO_save(std::string filename){
	IOFiles w(filename,true);
	for(unsigned int i(0);i<Nparticles_;i++){
		w<<px_[i]<<pv_[i]<<pbx_[i]<<pfbx_[i];
	}
}

void PSO::PSO_load(std::string filename){
	IOFiles r(filename,false);
	for(unsigned int i(0);i<Nparticles_;i++){
		r>>px_[i]>>pv_[i]>>pbx_[i]>>pfbx_[i];
	}
	for(unsigned int i(1);i<Nparticles_;i++){
		if(pfbx_[i] < pfbx_[bparticle_] ){
			bparticle_ = i;
		}
	}
}
/*}*/
