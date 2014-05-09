#include"PSO.hpp"
/*constructors and destructor*/
/*{*/
PSO::PSO(unsigned int Nbees, unsigned int Nfreedom, double cg, double cp, unsigned int maxiter):
	Nbees_(Nbees),
	Nfreedom_(Nfreedom),
	maxiter_(maxiter),
	pb_(new Vector<double>[Nbees]),
	pfb_(new double[Nbees]),
	bbee_(0),
	pv_(new Vector<double>[Nbees]),
	px_(new Vector<double>[Nbees]),
	min_(Nfreedom,-1.5),
	max_(Nfreedom,1.5),
	forget_(0),
	free_(new bool[Nbees]),
	cg_(cg),
	cp_(cp),
	chi_(-2.0/(2.0-(cp+cg)-sqrt((cp+cg)*(cp+cg)-4*(cp+cg)))),
	rnd_(1e4)
{
	if(cg_+cp_<4){
		std::cerr<<"PSO::PSO(Nparam,cg,cp,*f) : cg+cp<4 => chi=nan. Redefinition cg=cp=2.05"<<std::endl;
		cg_ = 2.05;
		cp_ = 2.05;
	}
	for(unsigned int i(0);i<Nbees_;i++){
		px_[i].set(Nfreedom_);
		pv_[i].set(Nfreedom_);
		pb_[i].set(Nfreedom_);
		free_[i] = true;
	}
}

void PSO::PSO_init(){
#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int i=0;i<Nbees_;i++){
#pragma omp critical
		{
			for(unsigned int j(0);j<Nfreedom_;j++){
				px_[i](j) = rnd_.get()*(max_(j)-min_(j))+min_(j);
				pv_[i](j) = rnd_.get()*(max_(j)-min_(j))+min_(j);
			}
		}
		pb_[i] = px_[i];
		pfb_[i] = run(px_[i]);
	}
	for(unsigned int i(1);i<Nbees_;i++){
		if(pfb_[i] < pfb_[bbee_] ){
			bbee_ = i;
		}
	}
}

void PSO::PSO_set_limit(unsigned int param, double min, double max){
	min_(param) = min;
	max_(param) = max;
}

PSO::~PSO(){
	if(pb_){delete[] pb_;}
	if(px_){delete[] px_;}
	if(pv_){delete[] pv_;}
	if(pfb_){delete[] pfb_;}
	if(free_){delete[] free_;}
}
/*}*/

/*core of the class*/
/*{*/
void PSO::move(unsigned int i){
	for(unsigned int j(0);j<Nfreedom_;j++){
		pv_[i](j) = chi_*(pv_[i](j) + cp_*rnd_.get()*(pb_[i](j)-px_[i](j)) + cg_*rnd_.get()*(pb_[bbee_](j)-px_[i](j)));
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

void PSO::evaluate(unsigned int i){
	double fx(run(px_[i]));
	if( fx < pfb_[i]){
		pfb_[i] = fx; 
		pb_[i] = px_[i]; 
	}
}

void PSO::launch(unsigned int i){
#pragma omp critical
	{
		move(i);
		if(forget_>maxiter_/Nbees_){
			std::cerr<<"PSO :: will forget "<< pb_[bbee_]<<std::endl;
			std::cerr<<"("<<forget_<<") with value  "<< pfb_[bbee_]<<std::endl;
			pb_[bbee_] = px_[bbee_];
			pfb_[bbee_] = 123456789.0;
			for(unsigned int j(0);j<Nbees_;j++){
				if(pfb_[j] < pfb_[bbee_] && j != bbee_){
					bbee_ = j;
				}
			}
			forget_ = 0;
		}
		forget_++;
	}
	evaluate(i);
	if(pfb_[i] < pfb_[bbee_] ){
		bbee_ = i;
		forget_ = 0;
	}
}

void PSO::PSO_run(bool synchro){
	std::cout<<"PSO::will start the run"<<std::endl;	
	if(synchro){
		for(unsigned int iter(0); iter<maxiter_; iter++){
#pragma omp parallel for
			for(unsigned int i=0;i<Nbees_;i++){
				launch(i);
			}
		}
	} else {
		//std::cerr<<"Nbees must be bigger than omp_get_num_threads"<<std::endl;
		unsigned int i(0);
		unsigned int ip(0);
#pragma omp parallel for schedule(dynamic,1) firstprivate(ip)
		for(unsigned int iter=0; iter<maxiter_; iter++){
#pragma omp critical
			{
				i++;
				ip = (i-1) % Nbees_;
				if(!free_[ip]){ iter--; }
			}
			if(free_[ip]){
				free_[ip]=false;
				launch(ip);
				free_[ip]=true;
			}
		}
	}
}
/*}*/

/*print, save and load*/
/*{*/
void PSO::PSO_print(){
	for(unsigned int i(0);i<Nbees_;i++){
		std::cout<<px_[i]<<pv_[i]<<pb_[i]<<pfb_[i]<<std::endl;
	}
}

void PSO::PSO_save(std::string filename){
	IOFiles w(filename,true);
	for(unsigned int i(0);i<Nbees_;i++){
		w<<px_[i]<<pv_[i]<<pb_[i]<<pfb_[i];
	}
}

void PSO::PSO_load(std::string filename){
	IOFiles r(filename,false);
	for(unsigned int i(0);i<Nbees_;i++){
		r>>px_[i]>>pv_[i]>>pb_[i]>>pfb_[i];
	}
	for(unsigned int i(1);i<Nbees_;i++){
		if(pfb_[i] < pfb_[bbee_] ){
			bbee_ = i;
		}
	}
}
/*}*/
