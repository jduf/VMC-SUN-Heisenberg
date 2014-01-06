#ifndef DEF_PSO
#define DEF_PSO

#include"Vector.hpp"
#include"Rand.hpp"
#include<omp.h>

class PSO {
	public:
		PSO(unsigned int Nbees, unsigned int Nfreedom, double cg, double cp);
		virtual ~PSO();

		void PSO_init(unsigned int bla);
		void PSO_set_limit(unsigned int param, double min, double max);
		void PSO_run(bool synchro=true);

	protected:
		virtual double run(Vector<double> const& x)=0;

		unsigned int Nbees_;
		unsigned int Nfreedom_;
		Vector<double>* pb_;//!< pb_[i] best position of the ith particle
		Vector<double> b_;//!< global best position
		double* pfb_; //!< pfb_[i] value at the best position of the ith particle
		double fb_;	//!< value at the global best position

	private:
		void launch(unsigned int i);
		void move(unsigned int i);
		void evaluate(unsigned int i);		

		unsigned int maxiter_;
		Vector<double>* pv_;//!< pv_[i] velocity of the ith particle
		Vector<double>* px_;//!< px_[i] position of the ith particle
		Vector<double> min_;//!< min_(c) minimum value of the cth coordinate
		Vector<double> max_; //!< max(c) minimum value of the cth coordinate
		bool* free_;//!< true if particle_[i] isn't running
		double cg_;
		double cp_;
		double chi_;
		Rand rnd_;
};
		
/*constructors and destructor*/
/*{*/
PSO::PSO(unsigned int Nbees, unsigned int Nfreedom, double cg, double cp):
	Nbees_(Nbees),
	Nfreedom_(Nfreedom),
	pb_(new Vector<double>[Nbees]),
	b_(Nfreedom),
	pfb_(new double[Nbees]),
	fb_(123456789.0),
	maxiter_(50),
	pv_(new Vector<double>[Nbees]),
	px_(new Vector<double>[Nbees]),
	min_(Nfreedom,-1.5),
	max_(Nfreedom,1.5),
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

void PSO::PSO_init(unsigned int bla){
	maxiter_ = bla;
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
	b_ = pb_[0];
	fb_= pfb_[0];
	for(unsigned int i(1);i<Nbees_;i++){
		if(pfb_[i] < fb_ ){
			b_ = pb_[i];
			fb_= pfb_[i];
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
		pv_[i](j) = chi_*(pv_[i](j) + cp_*rnd_.get()*(pb_[i](j)-px_[i](j)) + cg_*rnd_.get()*(b_(j)-px_[i](j)));
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
	}
	evaluate(i);
	if(pfb_[i] < fb_ ){
		fb_ = pfb_[i];
		b_ = pb_[i];
	}
}

void PSO::PSO_run(bool synchro){
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
				//if(ip<2){sleep(2);}
				//else{ sleep(5);}
				launch(ip);
				free_[ip]=true;
			}
		}
	}
}
/*}*/
#endif
