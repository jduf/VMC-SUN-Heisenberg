#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Lapack.hpp"
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
		/*! This method creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - sets N, m, n, sts_ and set tmp to the correct size 
		 * - allocates memory Ainv_
		 * - initialize the random number generator
		 * - creates an random initial state and computes its related matrices
		 */
		//}
		virtual unsigned int init(Container const& input, unsigned int thread);
		/*!exchanges two particles of different color */
		virtual void swap();
		/*!exchanges particle on site s1 with the one on site s2*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1);
		//{Description
		/*!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged, only one matrix is
		 *   modified, two of its columns are exchanged and therefore a minus
		 *   sign arises 
		 * - when two different colors are exchanged, computes the ratio using
		 *   the determinant lemma
		 */
		//}
		virtual Type ratio()=0;
		//{Description
		/*!Updates the state if the condition given by the System::ratio()
		 * method is accepted. The update consists of :
		 *
		 * - computes the Ainv_ matrices
		 * - updates the configuration : s
		 */
		//}
		
		virtual void update();
		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		void measure(double& E_config);
		virtual void print()=0;

		unsigned int N_;//!< number of different colors
		unsigned int n_;//!< number of lattice site
		unsigned int m_;//!< number of each color

		Matrix<unsigned int> s_;//!< on the i site : s(i,0)=color, s(i,1)=row

	protected:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		Rand* rnd;				//!< generator of random numbers 

		unsigned int color[2];	//!< colors of the exchanged sites
		unsigned int new_s[2];	//!< sites that are exchanged
		Matrix<unsigned int> sts_;//!< sts_(i,0) is a site that can be exchanged with sts_(i,1)
};

/*constructors and destructor*/
/*{*/
template<typename Type>
System<Type>::System():
	N_(0),
	n_(0),
	m_(0),
	rnd(NULL)
{ }

template<typename Type>
System<Type>::~System(){
	if(rnd){ delete rnd;}
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int System<Type>::init(Container const& input, unsigned int thread){
	N_ = input.get<unsigned int>("N");
	m_ = input.get<unsigned int>("m");
	n_ = input.get<unsigned int>("n");
	sts_ = input.get<Matrix<unsigned int> >("sts");

	rnd = new Rand(100,thread);
	s_.set(n_,2);

	return 1;
}

template<typename Type>
void System<Type>::update(){
	///*update the sites*/
	s_(new_s[0],0) = color[1];
	s_(new_s[1],0) = color[0];
}

template<typename Type>
void System<Type>::swap(){
	new_s[0] = rnd->get(n_);
	color[0] = s_(new_s[0],0);
	do {
		new_s[1] = rnd->get(n_);
		color[1] = s_(new_s[1],0);
	} while(color[0] == color[1]);
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1){
	/*if new_s[i] is not set, then System::jastrow() uses the wrong sites*/
	new_s[0] = s0;
	new_s[1] = s1;
	color[0] = s_(s0,0);
	color[1] = s_(s1,0);
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
void System<Type>::measure(double& E_config){
	E_config = 0.0;
	for(unsigned int i(0);i<sts_.row();i++){
		swap(sts_(i,0),sts_(i,1));
		/*!the minus sign is required because ∑<P_ij> is positive and this is
		 * what is actually computed*/
		E_config -= real(ratio());
	}
}
/*}*/
#endif
