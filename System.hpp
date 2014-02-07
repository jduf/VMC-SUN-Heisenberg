#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Rand.hpp"
#include "Container.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class System{
	public:
		/*!create a System without any parameters set*/
		System();

		/*!delete all the variables dynamically allocated*/
		virtual ~System();

		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - sets N, m, n, sts_ and set tmp to the correct size 
		 * - allocates memory Ainv_
		 * - initialize the random number generator
		 */ //}
		virtual unsigned int init(Container const& input, unsigned int const& thread);

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
		void measure(double& E_config);

		virtual void correlation(Matrix<unsigned int>* corr);
		virtual void print()=0;

		unsigned int n_;//!< sites' number
		unsigned int N_;//!< colors' number
		unsigned int m_;//!< particles per site' number
		unsigned int M_;//!< particles' number of each color

		Matrix<unsigned int> s_;//!< on the site i : s(i,0)=color, s(i,1)=row

	protected:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		Rand* rnd;				//!< generator of random numbers 

		unsigned int new_c[2];	//!< colors of the exchanged sites
		unsigned int new_s[2];	//!< sites that are exchanged
		unsigned int new_p[2];	//!< sites that are exchanged

	private:
		bool is_new_state_forbidden();

		Matrix<unsigned int> sts_;//!< sts_(i,0) is a site that can be exchanged with sts_(i,1)
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
System<Type>::System():
	n_(0),
	N_(0),
	m_(0),
	M_(0),
	rnd(NULL)
{ }

template<typename Type>
System<Type>::~System(){
	if(rnd){delete rnd;}
}

template<typename Type>
unsigned int System<Type>::init(Container const& input, unsigned int const& thread){
	n_ = input.get<unsigned int>("n");
	N_ = input.get<unsigned int>("N");
	m_ = input.get<unsigned int>("m");
	M_ = input.get<unsigned int>("M");

	sts_ = input.get<Matrix<unsigned int> >("sts");
	s_.set(n_,m_);

	rnd = new Rand(100,thread);
	return 1;
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
	new_s[0] = rnd->get(n_);
	new_p[0] = rnd->get(m_);
	new_c[0] = s_(new_s[0],new_p[0]);
	do {
		new_s[1] = rnd->get(n_);
		new_p[1] = rnd->get(m_);
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
void System<Type>::measure(double& E_config){
	E_config = 0.0;
	for(unsigned int i(0);i<sts_.row();i++){
		for(unsigned int p0(0); p0<m_; p0++){
			for(unsigned int p1(0); p1<m_; p1++){
				swap(sts_(i,0),sts_(i,1),p0,p1);
				if(!is_new_state_forbidden()){ E_config += real(ratio()); }
			}
		}
	}
}

template<typename Type>
void System<Type>::correlation(Matrix<unsigned int>* corr){if(corr){std::cerr<<"correlation not defined for System"<<std::endl;}}
/*}*/
#endif
