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
		/*!delete all the variables dynamically allocated with new*/
		~System();

		//{Description
		/*! This method creates the system in function of the input parameters.
		 *
		 * - sets N, n, sts, m and set tmp to the correct size 
		 * - allocates memory for s, A and Ainv
		 * - sets a different random number generator for each thread
		 * - creates an random initial state and computes its related matrices
		 */
		//}
		unsigned int init(Container const& input, unsigned int thread);
		/*!exchanges two particles of different color */
		void swap();
		/*!exchanges particle on site s1 with the one on site s2*/
		void swap(unsigned int const& s0, unsigned int const& s1);
		//{Description
		/*!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 *
		 * - when one matrix is modified, two of its columns are exchanged and
		 *   therefore a minus sign arises 
		 * - when two different color are exchanged, one computes the ratio
		 *   using the determinant lemma
		 */
		//}
		Type ratio();
		//{Description
		/*!Updates the state if the condition given by the System::ratio()
		 * method is acolorepted. The update consists of :
		 *
		 * - computes the new Ainv and A matrices
		 * - updates the configuration : s
		 */
		//}
		void update();
		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		void measure(double& E_config);
		void print();

		unsigned int N_;//!< number of different spin colors
		unsigned int n_;//!< number of lattice site
		unsigned int m_;//!< number of each color

	private:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		Rand* rnd;				//!< generator of random numbers 
		Matrix<Type> EVec;		//!< det(A) <=> <GS|a>
		Matrix<Type> *Ainv;		//!< inverse of A
		Matrix<Type> tmp;		//!< temporary matrix used during the update 
		Type w[2];				//!< det(W)= d = determinant ratios of <GS|a>/<GS|b> ; W=(w11,0,0,w22)
		Type d;					//!< Det(W)

		unsigned int color[2];
		unsigned int row[2];
		unsigned int new_ev[2];
		unsigned int new_s[2];
		Matrix<unsigned int> ev;
		Matrix<unsigned int> s;
		Matrix<unsigned int> sts;//!< sts(i,0) is a site that can be exchanged with sts(i,1)
};

/*constructors and destructor*/
/*{*/
template<typename Type>
System<Type>::System():
	N_(0),
	n_(0),
	m_(0),
	rnd(NULL),
	Ainv(NULL),
	d(0.0)
{ }

template<typename Type>
System<Type>::~System(){
	delete[] Ainv;
	delete rnd;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int System<Type>::init(Container const& input, unsigned int thread){
	N_ = input.get<unsigned int>("N");
	m_ = input.get<unsigned int>("m");
	n_ = input.get<unsigned int>("n");
	std::cout<<N_<<" "<<m_<<" "<<n_<<std::endl;
	sts = input.get<Matrix<unsigned int> >("sts");
	EVec= input.get<Matrix<Type> >("EVec");

	Ainv = new Matrix<Type>[N_];
	rnd = new Rand(100,thread);
	tmp.set(m_,m_);
	ev.set(N_,m_);
	s.set(n_,2);

	for(unsigned int i(0); i < N_; i++){
		Ainv[i].set(m_,m_);
	}

	unsigned int N_as(n_);
	unsigned int site(0);
	unsigned int* available(new unsigned int[n_]);
	Vector<int> P;
	unsigned int l(0);
	unsigned int TRY_MAX(100);
	double rcn(0.0);
	do {
		N_as = n_;
		for(unsigned int i(0); i < n_; i++){
			available[i] = i;
		}

		for(unsigned int spin(0); spin < N_; spin++){
			for(unsigned int i(0); i < m_; i++){
				site = rnd->get(N_as);
				ev(spin,i) = spin*n_+available[site];
				s(available[site],0) = spin;
				s(available[site],1) = i;
				for(unsigned int j(0); j < m_; j++){
					Ainv[spin](i,j) = EVec(ev(spin,i),j);
				}
				for(unsigned int j(site); j+1 < N_as; j++){
					available[j] = available[j+1];
				}
				N_as--;
			}
			Lapack<Type> inv(&Ainv[spin],false,'G');
			P = inv.is_singular(rcn);
			if(!P.ptr()){
				spin = N_;
			} else {
				inv.inv(P);
			}
		}
	} while (!P.ptr() && ++l<TRY_MAX);

	delete[] available;

	if(l==TRY_MAX){
		std::cerr<<"sorry, the thread will not be lunched because no initial state was found"<<std::endl;
		return 0;
	} else {
		std::cerr<<"yeah ! initial state found"<<std::endl;
		return 1;
	}
}

template<typename Type>
void System<Type>::swap(){
	new_s[0] = rnd->get(n_);
	color[0] = s(new_s[0],0);
	do {
		new_s[1] = rnd->get(n_);
		color[1] = s(new_s[1],0);
	} while(color[0]==color[1]);
	row[0] = s(new_s[0],1); 
	row[1] = s(new_s[1],1); 
	new_ev[0] = color[0]*n_ + new_s[1];
	new_ev[1] = color[1]*n_ + new_s[0];
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1){
	color[0] = s(s0,0);
	color[1] = s(s1,0);
	row[0] = s(s0,1);
	row[1] = s(s1,1); 
	new_ev[0] = color[0]*n_ + s1;
	new_ev[1] = color[1]*n_ + s0;
}

template<typename Type>
void System<Type>::update(){
	Type t_tmp;
	for(unsigned int c(0);c<2;c++){
		for(unsigned int j(0);j<m_;j++){
			if(row[c] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<m_;k++){
				t_tmp += EVec(new_ev[c],k)*Ainv[color[c]](k,j);
			}
			for(unsigned int i(0);i<m_;i++){
				tmp(i,j) = t_tmp*Ainv[color[c]](i,row[c])/w[c];
			}
		}
		Ainv[color[c]] -= tmp;
		ev(color[c],row[c]) = new_ev[c];
	}
	s(new_s[0],0) = color[1];
	s(new_s[0],1) = row[1];
	s(new_s[1],0) = color[0];
	s(new_s[1],1) = row[0];
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type System<Type>::ratio(){
	if(color[0]==color[1]){
		return -1.0;
	} else {
		for(unsigned int i(0);i<2;i++){
			w[i] = 0.0;
			for(unsigned int k(0);k<m_;k++){
				w[i] += EVec(new_ev[i],k)*Ainv[color[i]](k,row[i]);
			}
		}
		d=w[0]*w[1];
		if( std::abs(d) < 1e-10 ){ return 0.0; }
		else { return d; }
	}
}

template<typename Type>
void System<Type>::measure(double& E_config){
	E_config = 0.0;
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		/*the minus sign comes from <C|H|C'>*/
		E_config -= real(ratio());
	}
}

template<typename Type>
void System<Type>::print(){
	std::cout<<"=========================="<<std::endl;
	for(unsigned int spin(0);spin<N_;spin++){
		Matrix<Type> A(m_,m_);
		for(unsigned int i(0);i<m_;i++){
			for(unsigned int j(0);j<m_;j++){
				A(i,j) = EVec(ev(spin,i),j);
			}
		}
		std::cout<<(A*Ainv[spin]).diag().chop()<<std::endl;
	}
	std::cout<<"EVec"<<std::endl;
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0); i<m_; i++){
			std::cout<<ev(spin,i)<<" ";
		}
		std::cout<<std::endl;
	}
	for(unsigned int i(0); i < n_; i++){
		std::cout<<"("<<s(i,0)<<","<<s(i,1)<<") ";
		//std::cout<<s(i,0)<<" ";
	}
	std::cout<<std::endl;
	std::cout<<"=========================="<<std::endl;
}
/*}*/
#endif
