#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Lapack.hpp"
#include "Rand.hpp"

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
		 * - sets N_spin, N_site, sts, N_m and set tmp to the correct size 
		 * - allocates memory for s, A and Ainv
		 * - sets a different random number generator for each thread
		 * - creates an random initial state and computes its related matrices
		 */
		//}
		unsigned int init(unsigned int N_spin_, unsigned int N_m_, Matrix<unsigned int> const& sts_, Matrix<Type> const& EVec, unsigned int thread);
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
		void measure(double& E_config, Matrix<unsigned int>& lattice);
		void print();

		unsigned int N_spin;//!< number of different spin colors
		unsigned int N_site;//!< number of lattice site
		unsigned int N_m;	//!< number of each color

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
	N_spin(0),
	N_site(0),
	N_m(0),
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

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type System<Type>::ratio(){
	if(color[0]==color[1]){
		return -1.0;
	} else {
		for(unsigned int i(0);i<2;i++){
			w[i] = 0.0;
			for(unsigned int k(0);k<N_m;k++){
				w[i] += EVec(new_ev[i],k)*Ainv[color[i]](k,row[i]);
			}
		}
		d=w[0]*w[1];
		if( std::abs(d) < 1e-10 ){ return 0.0; }
		else { return d; }
	}
}

template<typename Type>
void System<Type>::measure(double& E_config, Matrix<unsigned int>& lattice){
	E_config = 0.0;
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		E_config -= real(ratio());
	}
	for(unsigned int i(0); i < N_site; i++){
		lattice(i,s(i,0)) += 1;
	}
}

template<typename Type>
void System<Type>::print(){
	std::cout<<"=========================="<<std::endl;
	for(unsigned int spin(0);spin<N_spin;spin++){
		Matrix<Type> A(N_m,N_m);
		for(unsigned int i(0);i<N_m;i++){
			for(unsigned int j(0);j<N_m;j++){
				A(i,j) = EVec(ev(spin,i),j);
			}
		}
		std::cout<<(A*Ainv[spin]).diag().transpose().chop()<<std::endl;
	}
	std::cout<<"EVec"<<std::endl;
	for(unsigned int spin(0);spin<N_spin;spin++){
		for(unsigned int i(0); i<N_m; i++){
			std::cout<<ev(spin,i)<<" ";
		}
		std::cout<<std::endl;
	}
	for(unsigned int i(0); i < N_site; i++){
		std::cout<<"("<<s(i,0)<<","<<s(i,1)<<") ";
		//std::cout<<s(i,0)<<" ";
	}
	std::cout<<std::endl;

	std::cout<<"=========================="<<std::endl;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int System<Type>::init(unsigned int N_spin_, unsigned int N_m_, Matrix<unsigned int> const& sts_, Matrix<Type> const& EVec_, unsigned int thread){
	N_spin = N_spin_;
	N_m = N_m_;
	N_site = N_spin*N_m;
	sts = sts_;
	EVec = EVec_;

	Ainv = new Matrix<Type>[N_spin];
	rnd = new Rand(100,thread);
	tmp.set(N_m,N_m);
	ev.set(N_spin,N_m);
	s.set(N_site,2);

	for(unsigned int i(0); i < N_spin; i++){
		Ainv[i].set(N_m,N_m);
	}

	unsigned int N_as(N_site);
	unsigned int site(0);
	unsigned int* available(new unsigned int[N_site]);
	Matrix<int> P;
	unsigned int l(0);
	unsigned int TRY_MAX(100);
	double rcn(0.0);
	do {
		N_as = N_site;
		for(unsigned int i(0); i < N_site; i++){
			available[i] = i;
		}

		for(unsigned int spin(0); spin < N_spin; spin++){
			for(unsigned int i(0); i < N_m; i++){
				site = rnd->get(N_as);
				ev(spin,i) = spin*N_site+available[site];
				s(available[site],0) = spin;
				s(available[site],1) = i;
				for(unsigned int j(0); j < N_m; j++){
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
				spin = N_spin;
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
		//print();
		return 1;
	}
}

template<typename Type>
void System<Type>::swap(){
	new_s[0] = rnd->get(N_site);
	color[0] = s(new_s[0],0);
	do {
		new_s[1] = rnd->get(N_site);
		color[1] = s(new_s[1],0);
	} while(color[0]==color[1]);
	row[0] = s(new_s[0],1); 
	row[1] = s(new_s[1],1); 
	new_ev[0] = color[0]*N_site + new_s[1];
	new_ev[1] = color[1]*N_site + new_s[0];
	//std::cout<<"si co row new_ev"<<std::endl;
	//std::cout<<new_s[0]<<" "<<color[0]<<" "<<row[0]<<" "<<new_ev[0]<<std::endl;
	//std::cout<<new_s[1]<<" "<<color[1]<<" "<<row[1]<<" "<<new_ev[1]<<std::endl;
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1){
	color[0] = s(s0,0);
	color[1] = s(s1,0);
	row[0] = s(s0,1);
	row[1] = s(s1,1); 
	new_ev[0] = color[0]*N_site + s1;
	new_ev[1] = color[1]*N_site + s0;
}

template<typename Type>
void System<Type>::update(){
	Type t_tmp;
	for(unsigned int m(0);m<2;m++){
		for(unsigned int j(0);j<N_m;j++){
			if(row[m] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<N_m;k++){
				t_tmp += EVec(new_ev[m],k)*Ainv[color[m]](k,j);
			}
			for(unsigned int i(0);i<N_m;i++){
				tmp(i,j) = t_tmp*Ainv[color[m]](i,row[m])/w[m];
			}
		}
		Ainv[color[m]] -= tmp;
		ev(color[m],row[m]) = new_ev[m];
	}
	s(new_s[0],0) = color[1];
	s(new_s[0],1) = row[1];
	s(new_s[1],0) = color[0];
	s(new_s[1],1) = row[0];

	//print();
}
/*}*/
#endif
