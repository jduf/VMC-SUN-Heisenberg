#ifndef DEF_MCSYSTEM
#define DEF_MCSYSTEM

#include "System.hpp"
#include "Rand.hpp"

/*!Abstract class that is used by MonteCarlo.hpp to sample the system.*/
template<typename Type>
class MCSystem: public virtual System{
	public:
		/*!Constructor*/
		MCSystem(System const& S, Rand& seed);
		/*!Destructor*/
		virtual ~MCSystem(){}

		/*!Exchanges two particles of different colors on random sites*/
		virtual void swap();
		/*!Exchanges the p0's particle of site s0 with the p1's of site s1*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);
		/*!Pure virtual method that computes the ratio between two states*/
		virtual Type ratio()=0;
		/*!Updates only s_*/
		virtual void update();

		/*!Sample the system for the new step*/
		void measure_new_step();
		/*!Add the sample to the statistic*/
		void add_sample();
		///*!Calls complete_analysis of the sampled datas*/
		void complete_analysis(double const& tol);
		
	protected:
		unsigned int new_c[2];//!< colors of the exchanged sites
		unsigned int new_s[2];//!< sites that are exchanged
		unsigned int new_p[2];//!< sites that are exchanged

		Matrix<unsigned int> s_;	//!< s(site,particle)=color
		Rand rnd_;					//!< generator of random numbers 

	private:
		/*!Forbid copy*/
		MCSystem(MCSystem<Type> const& mc);
		/*!Forbid assigment*/
		MCSystem& operator=(MCSystem<Type> const& mc);

		/*!Check only if the new state has not the same color on one site*/
		bool is_new_state_forbidden();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
MCSystem<Type>::MCSystem(System const& S, Rand& seed):
	System(S),
	s_(n_,m_),
	rnd_(1e3,seed)
{}
/*}*/

/*public method*/
/*{*/
template<typename Type>
void MCSystem<Type>::swap(){
	new_s[0] = rnd_.get(n_);
	new_p[0] = rnd_.get(m_);
	new_c[0] = s_(new_s[0],new_p[0]);
	do {
		new_s[1] = rnd_.get(n_);
		new_p[1] = rnd_.get(m_);
		new_c[1] = s_(new_s[1],new_p[1]);
	} while(is_new_state_forbidden() || new_c[0] == new_c[1]);
}

template<typename Type>
void MCSystem<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	new_s[0] = s0;
	new_p[0] = p0;
	new_c[0] = s_(s0,p0);
	new_s[1] = s1;
	new_p[1] = p1;
	new_c[1] = s_(s1,p1);
}

template<typename Type>
void MCSystem<Type>::update(){
	s_(new_s[0],new_p[0]) = new_c[1];
	s_(new_s[1],new_p[1]) = new_c[0];
}

template<typename Type>
void MCSystem<Type>::measure_new_step(){
	E_.set_x(0.0);
	double r;
	for(unsigned int i(0);i<links_.row();i++){
		corr_[i].set_x(0.0);
		for(unsigned int p0(0); p0<m_; p0++){
			for(unsigned int p1(0); p1<m_; p1++){
				swap(links_(i,0),links_(i,1),p0,p1);
				/*if the new state is forbidden, r=0 and therefore there is no
				 * need to complete the else condition*/
				if(!is_new_state_forbidden()){ 
					r = real(ratio());
					E_.add(r); 
					corr_[i].add(r);
				}
			}
		}
	}
	E_.divide(n_);

	if(lr_corr_.size()){
		for(unsigned int i(0);i<lr_corr_.size();i++){
			lr_corr_[i].set_x(0.0);
			/*might be a good idea to speed up this to take s over all n_*/
			for(unsigned int s(0);s<N_/m_;s++){
				for(unsigned int p0(0); p0<m_; p0++){
					for(unsigned int p1(0); p1<m_; p1++){
						swap(s,(i+s)%n_,p0,p1);
						if(!is_new_state_forbidden() && new_c[0] == new_c[1]){ lr_corr_[i].add(1); }
					}
				}
			}
		}
	}
}

template<typename Type>
void MCSystem<Type>::add_sample(){
	E_.add_sample();
	corr_.add_sample();
	lr_corr_.add_sample();
}

template<typename Type>
void MCSystem<Type>::complete_analysis(double const& tol){ 
	E_.complete_analysis(tol); 
	corr_.complete_analysis(tol); 
	lr_corr_.complete_analysis(tol); 
	for(unsigned int i(0);i<lr_corr_.size();i++){
		/*As the long range correlation are computed from N_/m_ different
		 * origins, one must compute the mean value*/
		lr_corr_[i].divide(N_/m_);
		lr_corr_[i].substract(1.0*m_*m_/N_);
	}
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
bool MCSystem<Type>::is_new_state_forbidden(){
	for(unsigned int i(0); i<m_; i++){
		if(s_(new_s[0],i) == new_c[1] && i != new_p[0]){ return true; }
		if(s_(new_s[1],i) == new_c[0] && i != new_p[1]){ return true; }
	}
	return false;
}
/*}*/
#endif
