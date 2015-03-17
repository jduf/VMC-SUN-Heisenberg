#include"PSO.hpp"

/*constructors and destructor*/
/*{*/
PSO::PSO(unsigned int Nparticles, unsigned int Nfreedom, double cg, double cp, unsigned int maxiter):
	Nparticles_(Nparticles),
	bparticle_(0),
	Nfreedom_(Nfreedom),
	maxiter_(maxiter),
	min_(Nfreedom,-2.5),
	max_(Nfreedom,2.5),
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
	/*can't use rnd_ inside a parallel region*/
	for(unsigned int i(0);i<Nparticles_;i++){
		for(unsigned int j(0);j<Nfreedom_;j++){
			px_[i](j) = rnd_.get()*(max_(j)-min_(j))+min_(j);
			pv_[i](j) = rnd_.get()*(max_(j)-min_(j))+min_(j);
		}
	}
#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int i=0;i<Nparticles_;i++){
		move_on_grid(i);
		pbx_[i] = px_[i];
		pfbx_[i] = f(px_[i]);
	}
	/*as bparticle_=0, start at i=1*/
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
	/*can't use rnd_ at the same time*/
#pragma omp critical(PSO_move_particle)
	{
		for(unsigned int j(0);j<Nfreedom_;j++){
			pv_[i](j) = chi_*(pv_[i](j) + cp_*rnd_.get()*(pbx_[i](j)-px_[i](j)) + cg_*rnd_.get()*(pbx_[bparticle_](j)-px_[i](j)));
			//if(v_[i] > (max_[i]-min_[i])/2.0){v_[i] = (max_[i]-min_[i])/4.0;}
			//if(v_[i] < (min_[i]-max_[i])/2.0){v_[i] = (min_[i]-max_[i])/4.0;}
			if( px_[i](j)+pv_[i](j) > max_(j)){ pv_[i](j) = log(1.0+rnd_.get()*(exp(max_(j)-px_[i](j))-1.0)); }
			if( px_[i](j)+pv_[i](j) < min_(j)){ pv_[i](j) =-log(1.0+rnd_.get()*(exp(px_[i](j)-min_(j))-1.0)); }
			px_[i](j) += pv_[i](j); 
		}
	}
}

void PSO::move_on_grid(unsigned int i){
	unsigned int n;
	double dx(0.05);
	for(unsigned int j(0);j<Nfreedom_;j++){
		n=0;
		if(std::abs(px_[i](j))<dx/2){ n=1; }
		if(std::abs(px_[i](j)-min_(j))<dx/2){ n=2; }
		if(std::abs(px_[i](j)-max_(j))<dx/2){ n=3; }
		switch(n){
			case 0:{ px_[i](j) = std::round(px_[i](j)/dx)*dx; }break;
			case 1:{ px_[i](j) = 0; }break;
			case 2:{ px_[i](j) = min_(j); }break;
			case 3:{ px_[i](j) = max_(j); }break;
		}
	}
#pragma omp critical
	{
		std::cout<<i<<":"<<n<<"->"<<px_[i]<<std::endl;
	}
}

void PSO::evaluate(unsigned int i){
	double fx(f(px_[i]));
	for(unsigned int i(0);i<Nparticles_;i++){
		if( are_equal(pbx_[i],px_[i]) ){ pfbx_[i] = fx; }
	}

	if( fx < pfbx_[i] ){
		pfbx_[i] = fx; 
		pbx_[i] = px_[i]; 
	}
	for(unsigned int j(0);j<Nparticles_;j++){
		if( are_equal(pbx_[j],px_[i]) ){
#pragma omp critical
			{
				std::cout<<"updated same position for i="<<i<<" j="<<j<<std::endl;
			}
			pfbx_[j] = fx; 
		}
	}
	if( pfbx_[i] < pfbx_[bparticle_] ){
#pragma omp critical
			{
				std::cout<<"best particle is now"<<i<<std::endl;
			}
		bparticle_ = i;
	}
}

void PSO::next_step(unsigned int i){
	move(i);
	move_on_grid(i);
	evaluate(i);
}

void PSO::PSO_run(){
	if(int(Nparticles_)<=omp_get_max_threads()){
		std::cout<<"PSO::run all particles in parallel"<<std::endl;
		for(unsigned int iter(0); iter<maxiter_; iter++){
#pragma omp parallel for
			for(unsigned int i=0;i<Nparticles_;i++){ next_step(i); }
		}
	} else {
		std::cout<<"PSO::run only "+tostring(omp_get_max_threads())+" particles in parallel"<<std::endl;
		unsigned int p(0);
#pragma omp parallel for schedule(dynamic,1) firstprivate(p)
		for(unsigned int i=0; i<maxiter_*Nparticles_; i++){
			p = (p+1) % Nparticles_;
			if(free_[p]){
				free_[p]=false;
				next_step(p);
				free_[p]=true;
			} else { i--; }
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
	std::cout<<"best at "<<pbx_[bparticle_]<<" with "<<pfbx_[bparticle_]<<std::endl;
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
