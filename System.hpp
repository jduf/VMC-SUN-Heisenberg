#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Matrice.hpp"
#include "Array2D.hpp"
#include "Lapack.hpp"
#include "Rand.hpp"
#include "Write.hpp"

#include <omp.h>
#include <cmath>

/*!Class that contains the information on the state
 * 
 *
 * 
*/
template<typename T>
class System{
	public:
		/*!create a System without any parameters set*/
		System();
		/*!delete all the variables dynamically allocated with new*/
		~System();

		//{Description
		/*! This method creates the system in function of the input parameters.
		 * The steps are :
		 *
		 * - sets N_spin, N_m, N_site, H and sts
		 * - allocates memory for A, Ainv and wis
		 * - sest a different random number generator for each thread
		 * - creates an random initial state and computes its related matrices
		 * - set w, cc, mc and tmp_m
		 */
		//}
		void init(unsigned int N_spin_, unsigned int N_m_, Matrice<double> const& H_, Array2D<unsigned int> const& sts_, Matrice<T> const& EVec, unsigned int thread);
		/*!exchanges two particles of different color */
		void swap();
		/*!exchanges particle on site s1 with the one on site s2*/
		void swap(unsigned int const& s1, unsigned int const& s2);
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
		T ratio();
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

		//Vecteur<double> bound;

	private:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assignment operator*/
		System& operator=(System const& S);

		Rand* rnd;			//!< generator of random numbers 
		Matrice<T> *A;      //!< det(A) <=> <GS|a>
		Matrice<T> *Ainv;   //!< inverse of A
		Matrice<T> tmp_m[2];//!< temporary matrices used during the update 
		T w[2];             //!< determinant ratios : <GS|a>/<GS|b>
		unsigned int a,b;	//!< two exchanged site by swap()
		unsigned int *wis;  //!< wis[i] = j : on ith site there is the j particle
		unsigned int mc[2]; //!< matrices (colors) that are modified 
		unsigned int cc[2]; //!< column's matrices (~band) that are exchanged 
		Matrice<double> H;	//!< Hamiltonian
		Array2D<unsigned int> sts; //!< sts(i,0) is a site that can be exchanged with sts(i,1)
};

/*constructors and destructor*/
/*{*/
template<typename T>
System<T>::System():
	N_spin(0),
	N_m(0),
	N_site(0),
	//bound(108,0),
	rnd(NULL),
	A(NULL),
	Ainv(NULL),
	a(0),
	b(0),
	wis(NULL),
	H(0,0),
	sts(0,0)
{ }

template<typename T>
System<T>::~System(){
	//for(unsigned int i(0);i<N_spin;i++){
		//Matrice<T> Ainv_check(A[i]);
		//Lapack<T> A_(Ainv_check.ptr(),Ainv_check.size(),'G');
		//A_.inv();
		//Matrice<T> check(Ainv_check*A[i]);
		//T t(0);
		//for(unsigned int j(0);j<check.size();j++){
			//t+= std::abs(check(j,j))-1.0;
		//}
		//T m(0);
		//for(unsigned int j(0);j<check.size();j++){
			//for(unsigned int k(0);k<check.size();k++){
				//m+= std::abs(Ainv[i](j,k)-Ainv_check(j,k));
			//}
		//}
		//std::cout<< i <<": trace of Ainv.A="<<t<<std::endl;
		//std::cout<< i <<": difference between Ainv and Ainv_check="<<m<<std::endl;
	//}

	delete[] A;
	delete[] Ainv;
	delete[] wis;
	delete rnd;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename T> 
T System<T>::ratio(){
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

template<typename T>
double System<T>::compute_energy(){
	double E_step(0.0);
	//double a(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		E_step += real(ratio() * H(sts(j,0),sts(j,1)));
		//a = -real(ratio());
		//E_step += a;
		//if(std::abs(H(sts(j,0),sts(j,1))+1)<1e-4){//if true -> th
			//bound(j) += a;
		//} else {
			//bound(j) += a;
		//}
	}
	return E_step;
}

template<typename T>
void System<T>::print(){
	for(unsigned int i(0); i<N_site; i++){
	std::cout<<wis[i]<<" ";
	}
	std::cout<<std::endl;
	//for(unsigned int i(0);i<2;i++){
		//std::cout<<mc[i]<<" "<< cc[i]<<" "<<w[i]<<std::endl;
	//}
	//for(unsigned int i(0);i<N_spin;i++){
		//std::cout<<std::endl;
		//std::cout<<Ainv[i] * A[i]<<std::endl;
	//}
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename T>
void System<T>::init(unsigned int N_spin_, unsigned int N_m_, Matrice<double> const& H_, Array2D<unsigned int> const& sts_, Matrice<T> const& EVec, unsigned int thread){
	N_spin = N_spin_;
	N_m = N_m_;
	N_site = N_spin*N_m;

	H = H_;
	sts = sts_;

	A = new Matrice<T>[N_spin];
	Ainv = new Matrice<T>[N_spin];
	wis = new unsigned int[N_site];
	rnd = new Rand(10,thread);

	unsigned int site(0);
	unsigned int* available_sites(new unsigned int[N_site]);
	unsigned int N_as(0);

	N_as = N_site;
	for(unsigned int i(0); i < N_as; i++){
		available_sites[i]  = i;	
	}

	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrice<T> (N_m);
		Ainv[i] = Matrice<T> (N_m);
		for(unsigned int j(0); j < N_m; j++){
			site = rnd->get(N_as);
			wis[available_sites[site]] = i*N_m+j; 	
			for(unsigned int k(0); k < N_m; k++){
				A[i](k,j) = EVec(available_sites[site],k);
			}
			for(unsigned int k(site); k < N_as-1; k++){
				available_sites[k] = available_sites[k+1];
			}
			N_as--;
		}
		Ainv[i] = A[i];
		Lapack<T> A_(Ainv[i].ptr(),Ainv[i].size(),'G');
		A_.inv();
		//Matrice<T> check(Ainv[i]*A[i]);
		//T t(0);
		//for(unsigned int j(0);j<check.size();j++){
			//t+= std::abs(check(j,j))-1.0;
		//}
		//std::cout<< i <<": trace of Ainv.A="<<t<<std::endl;
	}

	delete[] available_sites;

	for(unsigned int i(0);i<2;i++){
		w[i] = 0;
		cc[i] = 0;
		mc[i] = 0;
		tmp_m[i] = Matrice<T>(N_m);
	}
}

template<typename T>
void System<T>::swap() {
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

template<typename T>
void System<T>::swap(unsigned int const& s1, unsigned int const& s2) {
	mc[0] = wis[s1] / N_m; //gives the color of particle on site s1
	mc[1] = wis[s2] / N_m; //gives the color of particle on site s2
	cc[0] = wis[s1] % N_m; //gives the band of particle on site s1
	cc[1] = wis[s2] % N_m; //gives the band of particle on site s2
}

template<typename T>
void System<T>::update(){
	unsigned int s_tmp(wis[b]);
	wis[b] = wis[a];
	wis[a] = s_tmp;

	// there is a way to avoid this loop and its useless copy if one keeps track
	// only of the columns that are echanged...
	// exchange the two columns
	T tmp(0.0); 
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

/*double real(T)*/
/*{*/
template<typename T>
double real(T const& x);

template<>
inline double real(double const& x){
	return x;
}

template<>
inline double real(std::complex<double> const& x){
	return std::real(x);
}
/*}*/
#endif
