#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "CreateSystem.hpp"
#include "Rand.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class System{
	public:
		/*!Create a System and extract the parameters from CreateSystem*/
		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - sets N, m, n, links_ and set tmp to the correct size 
		 * - allocates memory Ainv_
		 * - initialize the random number generator
		 */ //}
		System(CreateSystem* CS, unsigned int const& thread);
		/*!Delete all the variables dynamically allocated*/
		virtual ~System();

		/*!Exchanges two particles of different color */
		virtual void swap();
		/*!Exchanges particle on site s1 with the one on site s2*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);
		/*!Virtual method that is called by MonteCarlo */
		virtual Type ratio()=0;
		/*!Updates only s_*/
		virtual void update();

		/*!Returns the status*/
		unsigned int get_status() const {return status_;}

		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		void measure(unsigned int type);	
		double get_energy_step() const {return E_step_;};
		Vector<double> get_corr_step() const {return corr_step_;};
		Vector<double> get_long_range_corr_step() const {return long_range_corr_step_;};

		/*!Pure virtual function that provides a way to check the System */
		virtual void print()=0;

	protected:
		unsigned int new_c[2];	//!< colors of the exchanged sites
		unsigned int new_s[2];	//!< sites that are exchanged
		unsigned int new_p[2];	//!< sites that are exchanged

		unsigned int status_;//!< Degenerate :0 No initial state:1

		unsigned int const N_;//!< colors' number
		unsigned int const n_;//!< sites' number
		unsigned int const m_;//!< particles per site' number
		unsigned int const M_;//!< particles' number of each color

		Matrix<unsigned int> s_;//!< on the site i : s(i,0)=color, s(i,1)=row
		Matrix<unsigned int> const links_;	//!< list of links

		Rand* rnd_;			//!< generator of random numbers 

		double E_step_;
		Vector<double> corr_step_;	//!< Correlation for each link 
		Vector<double> long_range_corr_step_;//!< Correlation for each link 

	private:
		/*!Check only if the new state has not the same color on one site*/
		bool is_new_state_forbidden();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
System<Type>::System(CreateSystem* CS, unsigned int const& thread):
	status_(1),
	N_(CS->get_N()),
	n_(CS->get_n()),
	m_(CS->get_m()),
	M_((m_*n_)/N_),
	s_(n_,m_),
	links_(CS->get_links()),
	rnd_(new Rand(100,thread)),
	E_step_(0.0),
	corr_step_(CS->get_num_links(),0),
	long_range_corr_step_(CS->get_n()/3-1,0)
{}

template<typename Type>
System<Type>::~System(){
	if(rnd_){delete rnd_;}
}
/*}*/

/*public method*/
/*{*/
template<typename Type>
void System<Type>::swap(){
	new_s[0] = rnd_->get(n_);
	new_p[0] = rnd_->get(m_);
	new_c[0] = s_(new_s[0],new_p[0]);
	do {
		new_s[1] = rnd_->get(n_);
		new_p[1] = rnd_->get(m_);
		new_c[1] = s_(new_s[1],new_p[1]);
	} while(is_new_state_forbidden() || new_c[0] == new_c[1]);
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	new_s[0] = s0;
	new_p[0] = p0;
	new_c[0] = s_(s0,p0);
	new_s[1] = s1;
	new_p[1] = p1;
	new_c[1] = s_(s1,p1);
}

template<typename Type>
void System<Type>::update(){
	/*update the sites*/
	s_(new_s[0],new_p[0]) = new_c[1];
	s_(new_s[1],new_p[1]) = new_c[0];
}

template<typename Type>
void System<Type>::measure(unsigned int type){
	E_step_ = 0.0;
	double r;
	for(unsigned int i(0);i<links_.row();i++){
		corr_step_(i) = 0.0;
		for(unsigned int p0(0); p0<m_; p0++){
			for(unsigned int p1(0); p1<m_; p1++){
				swap(links_(i,0),links_(i,1),p0,p1);
				if(!is_new_state_forbidden()){ 
					r = real(ratio());
					E_step_ += r; 
					corr_step_(i) += r;
				}
			}
		}
	}
	E_step_ /= n_;
	if(type==2){
		unsigned int x0(n_/3);
		for(unsigned int i(0);i<long_range_corr_step_.size();i++){
			long_range_corr_step_(i) = 0.0;
			for(unsigned int p0(0); p0<m_; p0++){
				for(unsigned int p1(0); p1<m_; p1++){
					swap(x0,x0+i+1,p0,p1);
					if(!is_new_state_forbidden()){ 
						long_range_corr_step_(i) += real(ratio());;
					}
				}
			}
		}
	}
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
bool System<Type>::is_new_state_forbidden(){
	for(unsigned int i(0); i<this->m_; i++){
		if(s_(new_s[0],i) == new_c[1] && i != new_p[0]){ return true; }
		if(s_(new_s[1],i) == new_c[0] && i != new_p[1]){ return true; }
	}
	return false;
}
/*}*/
#endif
