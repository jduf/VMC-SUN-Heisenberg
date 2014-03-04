#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "CreateSystem.hpp"
#include "Rand.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class System{
	public:
		/*!create a System without any parameters set*/
		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - sets N, m, n, sts_ and set tmp to the correct size 
		 * - allocates memory Ainv_
		 * - initialize the random number generator
		 */ //}
		System(CreateSystem* CS, unsigned int const& thread);

		/*!delete all the variables dynamically allocated*/
		virtual ~System();

		/*!Exchanges two particles of different color */
		virtual void swap();

		/*!Exchanges particle on site s1 with the one on site s2*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);

		/*!Virtual method that is called by MonteCarlo */
		virtual Type ratio()=0;

		//{Description
		/*!Updates the configuration s_ if the condition given by the
		 * System::ratio() method is accepted.
		 */ //}
		virtual void update();

		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		void measure(double& E_step, Vector<double>& corr_step);

		virtual void print()=0;

		unsigned int get_status() const {return status_;}
		unsigned int get_n() const {return n_;}

	protected:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		unsigned int new_c[2];	//!< colors of the exchanged sites
		unsigned int new_s[2];	//!< sites that are exchanged
		unsigned int new_p[2];	//!< sites that are exchanged

		unsigned int const N_;//!< colors' number
		unsigned int const n_;//!< sites' number
		unsigned int const m_;//!< particles per site' number
		unsigned int const M_;//!< particles' number of each color

		unsigned int status_;

		Matrix<unsigned int> s_;//!< on the site i : s(i,0)=color, s(i,1)=row
		Matrix<unsigned int> sts_;//!< sts_(i,0) is a site that can be exchanged with sts_(i,1)

		Rand* rnd_;				//!< generator of random numbers 

	private:
		bool is_new_state_forbidden();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
System<Type>::System(CreateSystem* CS, unsigned int const& thread):
	N_(CS->get_N()),
	n_(CS->get_n()),
	m_(CS->get_m()),
	M_((m_*n_)/N_),
	status_(1),
	s_(n_,m_),
	sts_(CS->get_sts()),
	rnd_(new Rand(100,thread))
{}

template<typename Type>
System<Type>::~System(){
	if(rnd_){delete rnd_;}

}

/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void System<Type>::update(){
	///*update the sites*/
	s_(new_s[0],new_p[0]) = new_c[1];
	s_(new_s[1],new_p[1]) = new_c[0];
}

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
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
bool System<Type>::is_new_state_forbidden(){
	for(unsigned int i(0); i<this->m_; i++){
		if(s_(new_s[0],i) == new_c[1] && i != new_p[0]){ return true; }
		if(s_(new_s[1],i) == new_c[0] && i != new_p[1]){ return true; }
	}
	return false;
}

template<typename Type>
void System<Type>::measure(double& E_step, Vector<double>& corr_step){
	E_step = 0.0;
	double r;
	for(unsigned int i(0);i<sts_.row();i++){
		corr_step(i) = 0.0;
		for(unsigned int p0(0); p0<m_; p0++){
			for(unsigned int p1(0); p1<m_; p1++){
				swap(sts_(i,0),sts_(i,1),p0,p1);
				if(!is_new_state_forbidden()){ 
					r = real(ratio());
					E_step += r; 
					corr_step(i) += r;
				}
			}
		}
	}
}
/*}*/
#endif
