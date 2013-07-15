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
		 * - sets N_spin, N_m, N_site, H and sts
		 * - allocates memory for A, Ainv and s
		 * - sets a different random number generator for each thread
		 * - creates an random initial state and computes its related matrices
		 * - set cc, sc, w and tmp
		 */
		//}
		unsigned int init(unsigned int N_spin_, unsigned int N_m_, Matrix<unsigned int> const& sts_, Matrix<Type> const& EVec, unsigned int thread);
		/*!exchanges two particles of different color */
		void swap();
		/*!exchanges particle on site s1 with the one on site s2*/
		void swap(unsigned int const& s0, unsigned int const& s1);
		//{Description
		/*!Computes the ratio of the two determinants related to the current and next
		 * configuration
		 *
		 * - when one matrix is modified, two of its columns are exchanged and
		 *   therefore a minus sign arises 
		 * - when two matrices are modified, one computes the ratio using the
		 *   determinant lemma
		 */
		//}
		Type ratio();
		//{Description
		/*!Updates the state if the condition given by the System::ratio()
		 * method is accepted. The update consists of :
		 *
		 * - computes the new matrices
		 * - computes the new inverse matrices with the Sherman-Morisson formula
		 * - updates the configuration (s)
		 */
		//}
		void update();
		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		double compute_energy();

		void print();

		unsigned int N_spin;//!< number of different spin colors
		unsigned int N_site;//!< number of lattice site
		unsigned int N_m;	//!< number of each color

	private:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		Rand* rnd;			//!< generator of random numbers 
		Matrix<Type> *A;    //!< det(A) <=> <GS|a>
		Matrix<Type> *Ainv;   //!< inverse of A
		Matrix<Type> tmp;     //!< temporary matrix used during the update 
		Type w[2];             //!< determinant ratios : <GS|a>/<GS|b>
		Type d;           //!< Det(W)
		unsigned int *s;    //!< s[i] = j : on ith site there is the j particle
		unsigned int new_s[2];//!< two sites exchanged by swap()
		unsigned int sc[2]; //!< sites that are exchanged
		unsigned int cc[2]; //!< colors that are exchanged 
		Matrix<unsigned int> sts; //!< sts(i,0) is a site that can be exchanged with sts(i,1)
};

/*double real(T)*/
/*{*/
template<typename Type>
double real(Type const& x);

template<>
inline double real(double const& x){
	return x;
}

template<>
inline double real(std::complex<double> const& x){
	return std::real(x);
}
/*}*/

/*constructors and destructor*/
/*{*/
template<typename Type>
System<Type>::System():
	N_spin(0),
	N_site(0),
	N_m(0),
	rnd(NULL),
	Ainv(NULL),
	d(0.0),
	s(NULL)
{ }

template<typename Type>
System<Type>::~System(){
	delete[] A;
	delete[] Ainv;
	delete[] s;
	delete rnd;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type System<Type>::ratio(){
	if(cc[0]==cc[1]){
		return -1.0;
	} else {
		w[0] = 0.0;
		w[1] = 0.0;
		for(unsigned int k=0;k<N_m;k++){
			w[0] += A[cc[1]](sc[1],k)*Ainv[cc[0]](k,sc[0]);
			w[1] += A[cc[0]](sc[0],k)*Ainv[cc[1]](k,sc[1]);
		}
		d=w[0]*w[1];
		if( std::abs(d) < 1e-10 ){ return 0.0; }
		else { return d; }
	}
}

template<typename Type>
double System<Type>::compute_energy(){
	double E_step(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		E_step -= real(ratio());
	}
	return E_step;
}

template<typename Type>
void System<Type>::print(){
	//for(unsigned int i(0); i<N_site; i++){
	//std::cout<<s[i]<<" ";
	//}
	//std::cout<<std::endl;
	//for(unsigned int i(0); i < N_spin; i++){
	//std::cout<<(Ainv[i]*A[i]).diag()<<std::endl;;
	//}
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int System<Type>::init(unsigned int N_spin_, unsigned int N_m_, Matrix<unsigned int> const& sts_, Matrix<Type> const& EVec, unsigned int thread){
	N_spin = N_spin_;
	N_m = N_m_;
	N_site = N_spin*N_m;
	sts = sts_;

	A = new Matrix<Type>[N_spin];
	Ainv = new Matrix<Type>[N_spin];
	tmp.set(N_m,N_m);
	s = new unsigned int[N_site];
	rnd = new Rand(10,thread);

	for(unsigned int i(0); i < N_spin; i++){
		A[i].set(N_m,N_m);
		Ainv[i].set(N_m,N_m);
	}

	unsigned int N_as(N_site);
	unsigned int spin(0);
	unsigned int* available_s(new unsigned int[N_site]);
	Matrix<int> P;
	unsigned int l(0);
	unsigned int TRY_MAX(100);
	double rcn(0.0);
	std::cerr<<"System : init : get rid of A[i] and save only Ainv[i] & EVec => modify swap and update"<<std::endl;
	do {
		N_as = N_site;
		for(unsigned int i(0); i < N_site; i++){
			available_s[i] = i;
		}

		for(unsigned int i(0); i < N_spin; i++){
			for(unsigned int j(0); j < N_m; j++){
				spin = rnd->get(N_as);
				s[available_s[spin]] = i*N_m+j;
				for(unsigned int k(0); k < N_m; k++){
					A[i](j,k) = EVec(available_s[spin],k);
				}
				for(unsigned int k(spin); k < N_as-1; k++){
					available_s[k] = available_s[k+1];
				}
				N_as--;
			}
			Lapack<Type> inv(&A[i],true,'G');
			P = inv.is_singular(rcn);
			if(!P.ptr()){
				i = N_spin;
			} else {
				inv.inv(P);
				Ainv[i] = inv.get_mat();
			}
		}
	} while (!P.ptr() && ++l<TRY_MAX);

	delete[] available_s;

	if(l==TRY_MAX){
		std::cerr<<"sorry, the thread will not be lunched because no initial state was found"<<std::endl;
		return 0;
	} else {
		std::cerr<<"yeah ! initial state found"<<std::endl;
		print();
		return 1;
	}
}

template<typename Type>
void System<Type>::swap(){
	new_s[0] = rnd->get(N_site);
	new_s[1] = rnd->get(N_site);
	cc[0] = s[new_s[0]] / N_m; //gives the color of particle on site a
	cc[1] = s[new_s[1]] / N_m; //gives the color of particle on site b
	while(cc[0]==cc[1]){
		new_s[1] = rnd->get(N_site);
		cc[1] = s[new_s[1]] / N_m; //gives the color of particle on site b
	}
	sc[0] = s[new_s[0]] % N_m; //gives the band of particle on site a
	sc[1] = s[new_s[1]] % N_m; //gives the band of particle on site b
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1) {
	cc[0] = s[s0] / N_m; //gives the color of particle on site s1
	cc[1] = s[s1] / N_m; //gives the color of particle on site s2
	sc[0] = s[s0] % N_m; //gives the band of particle on site s1
	sc[1] = s[s1] % N_m; //gives the band of particle on site s2
}

template<typename Type>
void System<Type>::update(){
	unsigned int i_tmp(s[new_s[1]]);
	s[new_s[1]] = s[new_s[0]];
	s[new_s[0]] = i_tmp;

	// there is a way to avoid this loop and its useless copy if one keeps track
	// only of the columns that are echanged...
	// exchange the two columns
	Type t_tmp(0.0); 
	for(unsigned int j(0); j<N_m; j++){
		t_tmp = A[cc[0]](sc[0],j);
		A[cc[0]](sc[0],j) = A[cc[1]](sc[1],j);
		A[cc[1]](sc[1],j) = t_tmp;
	}

	//compute Ainv
	for(unsigned int m=0;m<2;m++){
		for(unsigned int j(0);j<N_m;j++){
			if(sc[m] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<N_m;k++){
				t_tmp += A[cc[m]](sc[m],k)*Ainv[cc[m]](k,j);
			}
			for(unsigned int i(0);i<N_m;i++){
				tmp(i,j) = t_tmp*Ainv[cc[m]](i,sc[m])/w[m];
			}
		}
		Ainv[cc[m]] -= tmp;
	}
}
/*}*/
#endif
