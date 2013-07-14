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
		 * - allocates memory for A, Ainv and wis
		 * - sets a different random number generator for each thread
		 * - creates an random initial state and computes its related matrices
		 * - set w, cc, mc and tmp_m
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
		 * - updates the configuration (wis)
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
		unsigned int N_m;	//!< number of each color
		unsigned int N_site;//!< number of lattice site

	private:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		Rand* rnd;			//!< generator of random numbers 
		Matrix<Type> *A;      //!< det(A) <=> <GS|a>
		Matrix<Type> *Ainv;   //!< inverse of A
		Matrix<Type> tmp_m[2];//!< temporary matrices used during the update 
		Type w[2];             //!< determinant ratios : <GS|a>/<GS|b>
		unsigned int a,b;	//!< two exchanged site by swap()
		unsigned int *wis;  //!< wis[i] = j : on ith site there is the j particle
		unsigned int mc[2]; //!< matrices (colors) that are modified 
		unsigned int cc[2]; //!< column's matrices (~band) that are exchanged 
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
	N_m(0),
	N_site(0),
	rnd(NULL),
	A(NULL),
	Ainv(NULL),
	a(0),
	b(0),
	wis(NULL),
	sts(0,0)
{ }

template<typename Type>
System<Type>::~System(){
	delete[] A;
	delete[] Ainv;
	delete[] wis;
	delete rnd;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type System<Type>::ratio(){
	if(mc[0] == mc[1]){
		return -1.0;
	} else {
		w[0] = 0.0;
		w[1] = 0.0;
		for(unsigned int i=0;i<N_m;i++){
			w[0] += Ainv[mc[0]](cc[0],i)*A[mc[1]](i,cc[1]);
			w[1] += Ainv[mc[1]](cc[1],i)*A[mc[0]](i,cc[0]);
		}
		if( std::abs(w[0]*w[1]) < 1e-10 ){ return 0.0; }
		else { return w[0]*w[1]; }
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
	//std::cout<<wis[i]<<" ";
	//}
	//std::cout<<std::endl;
	for(unsigned int i(0); i < N_spin; i++){
		std::cout<<(Ainv[i]*A[i]).diag()<<std::endl;;
	}
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int System<Type>::init(unsigned int N_spin_, unsigned int N_m_, Matrix<unsigned int> const& sts_, Matrix<Type> const& EVec, unsigned int thread){
	N_m = N_m_;
	N_spin = N_spin_;
	N_site = N_spin*N_m;
	sts = sts_;

	A = new Matrix<Type>[N_spin];
	Ainv = new Matrix<Type>[N_spin];
	wis = new unsigned int[N_site];
	rnd = new Rand(10,thread);

	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrix<Type> (N_m,N_m);
		Ainv[i] = Matrix<Type> (N_m,N_m);
	}

	unsigned int N_available(N_site);
	unsigned int spin(0);
	unsigned int* available_spin(new unsigned int[N_site]);
	Matrix<int> P;
	unsigned int l(0);
	unsigned int TRY_MAX(100);
	double rcn(0.0);
	do {
		N_available = N_site;
		for(unsigned int i(0); i < N_site; i++){
			available_spin[i]  = i%N_spin;
		}

		for(unsigned int i(0); i < N_site; i++){
			spin = rnd->get(N_available);
			wis[i] = i+available_spin[spin]; 	
			std::cout<<"test "<<i<<" "<<available_spin[spin]<<" "<<wis[i]<<std::endl;
			for(unsigned int k(0); k < N_m; k++){
				A[available_spin[spin]](k,i%N_m) = EVec(wis[i],k);
			}
			for(unsigned int k(spin); k < N_available-1; k++){
				available_spin[k] = available_spin[k+1];
			}
			N_available--;
		}

		for(unsigned int i(0);i<N_spin;i++){
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

delete[] available_spin;

for(unsigned int i(0);i<2;i++){
	w[i] = 0;
	cc[i] = 0;
	mc[i] = 0;
	tmp_m[i] = Matrix<Type>(N_m,N_m);
}
if(l==TRY_MAX){
	std::cerr<<"sorry, the thread will not be lunched because no initial state was found"<<std::endl;
	return 0;
} else {
	std::cerr<<"yeah ! initial state found, rcn="<<rcn<<std::endl;
	print();
	return 1;
}
}

template<typename Type>
void System<Type>::swap() {
	a = rnd->get(N_site);
	b = rnd->get(N_site);
	mc[0] = wis[a] / N_m; //gives the color of particle on site a
	mc[1] = wis[b] / N_m; //gives the color of particle on site b
	while(mc[0]==mc[1]){
		b = rnd->get(N_site);
		mc[1] = wis[b] / N_m; //gives the color of particle on site b
	}
	cc[0] = wis[a] % N_m; //gives the band of particle on site a
	cc[1] = wis[b] % N_m; //gives the band of particle on site b
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1) {
	mc[0] = wis[s0] / N_m; //gives the color of particle on site s1
	mc[1] = wis[s1] / N_m; //gives the color of particle on site s2
	cc[0] = wis[s0] % N_m; //gives the band of particle on site s1
	cc[1] = wis[s1] % N_m; //gives the band of particle on site s2
}

template<typename Type>
void System<Type>::update(){
	unsigned int s_tmp(wis[b]);
	wis[b] = wis[a];
	wis[a] = s_tmp;

	// there is a way to avoid this loop and its useless copy if one keeps track
	// only of the columns that are echanged...
	// exchange the two columns
	Type tmp(0.0); 
	for(unsigned int i(0); i<N_m; i++){
		tmp = A[mc[0]](i,cc[0]);
		A[mc[0]](i,cc[0]) = A[mc[1]](i,cc[1]);
		A[mc[1]](i,cc[1]) = tmp;
	}

	//compute Ainv
	for(unsigned int m=0;m<2;m++){
		for(unsigned int i(0);i<N_m;i++){
			if(cc[m] == i){ tmp = -1.0; }
			else { tmp = 0.0; }
			for(unsigned int j(0);j<N_m;j++){
				tmp += Ainv[mc[m]](i,j)*A[mc[m]](j,cc[m]);
			}
			for(unsigned int j(0);j<N_m;j++){
				tmp_m[m](i,j) = tmp*Ainv[mc[m]](cc[m],j)/w[m];
			}
		}
		Ainv[mc[m]] -= tmp_m[m];
	}
}
/*}*/
#endif
