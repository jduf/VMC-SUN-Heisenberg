#ifndef DEF_MCSYSTEM
#define DEF_MCSYSTEM

#include "CreateSystem.hpp"
#include "Rand.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class MCSystem: public virtual System{
	public:
		MCSystem(System const& S, unsigned int const& type);
		virtual ~MCSystem();

		/*!Exchanges two particles of different colors*/
		virtual void swap();
		/*!Exchanges particle on site s1 with the one on site s2*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);
		/*!Virtual method that is called by MonteCarlo */
		virtual Type ratio()=0;
		/*!Updates only s_*/
		virtual void update();

		/*!Returns the status*/
		bool found_initial_state() const {return found_initial_state_;}

		/*{Description
		 * !Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation 
		}*/
		void measure_new_step();	
		void add_sample();
		bool is_converged(double const& tol);
		void complete_analysis(double const& tol);
		
		void set();
		void init(unsigned int const& thread);

	protected:
		unsigned int new_c[2];//!< colors of the exchanged sites
		unsigned int new_s[2];//!< sites that are exchanged
		unsigned int new_p[2];//!< sites that are exchanged

		unsigned int type_;

		bool found_initial_state_;	

		Rand* rnd_;	//!< generator of random numbers 

		Matrix<unsigned int> s_;//!< on the site i : s(i,0)=color, s(i,1)=row

	private:
		MCSystem(MCSystem<Type> const&);
		MCSystem& operator=(MCSystem<Type> const&);
		virtual void init() = 0;

		/*!Check only if the new state has not the same color on one site*/
		bool is_new_state_forbidden();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
MCSystem<Type>::MCSystem(System const& S, unsigned int const& type):
	System(S),
	type_(type),
	found_initial_state_(false),
	rnd_(NULL)
{}

template<typename Type>
MCSystem<Type>::~MCSystem(){
	if(rnd_){delete rnd_;}
}

template<typename Type>
void MCSystem<Type>::init(unsigned int const& thread){
	s_.set(this->n_,this->M_);
	rnd_ = new Rand(100,thread);
	set(); 
	init();
}
/*}*/

/*public method*/
/*{*/
template<typename Type>
void MCSystem<Type>::swap(){
	new_s[0] = rnd_->get(this->n_);
	new_p[0] = rnd_->get(this->m_);
	new_c[0] = s_(new_s[0],new_p[0]);
	do {
		new_s[1] = rnd_->get(this->n_);
		new_p[1] = rnd_->get(this->m_);
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
	/*update the sites*/
	s_(new_s[0],new_p[0]) = new_c[1];
	s_(new_s[1],new_p[1]) = new_c[0];
}

template<typename Type>
void MCSystem<Type>::measure_new_step(){
	E_ = 0.0;
	double r;
	for(unsigned int i(0);i<this->links_.row();i++){
		corr_[i] = 0.0;
		for(unsigned int p0(0); p0<this->m_; p0++){
			for(unsigned int p1(0); p1<this->m_; p1++){
				swap(links_(i,0),links_(i,1),p0,p1);
				if(!is_new_state_forbidden()){ 
					r = real(ratio());
					E_ += r; 
					corr_[i] += r;
				}
			}
		}
	}
	E_ /= this->n_;
	if(long_range_corr_.size()!=0){
		unsigned int x0(this->n_/3);
		for(unsigned int i(0);i<this->n_/3;i++){
			long_range_corr_[i] = 0.0;
			for(unsigned int p0(0); p0<this->m_; p0++){
				for(unsigned int p1(0); p1<this->m_; p1++){
					swap(x0,x0+i+1,p0,p1);
					if(!is_new_state_forbidden() && new_c[0]!=new_c[1]){ 
						long_range_corr_[i] += real(ratio());;
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
	long_range_corr_.add_sample();
}

template<typename Type>
bool MCSystem<Type>::is_converged(double const& tol){ 
	//corr_.compute_convergence(tol); 
	//long_range_corr_.compute_convergence(tol); 
	E_.compute_convergence(tol); 
	return E_.get_conv();
}

template<typename Type>
void MCSystem<Type>::complete_analysis(double const& tol){ 
	E_.complete_analysis(tol); 
	corr_.complete_analysis(tol); 
	long_range_corr_.complete_analysis(tol); 
}

template<typename Type>
void MCSystem<Type>::set(){
	E_.set(50,5,false);
	corr_.set(this->n_,50,5,false);
	if(type_ == 2){ long_range_corr_.set(this->n_/3,50,5,false); }
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
bool MCSystem<Type>::is_new_state_forbidden(){
	for(unsigned int i(0); i<this->m_; i++){
		if(s_(new_s[0],i) == new_c[1] && i != new_p[0]){ return true; }
		if(s_(new_s[1],i) == new_c[0] && i != new_p[1]){ return true; }
	}
	return false;
}
/*}*/
#endif

