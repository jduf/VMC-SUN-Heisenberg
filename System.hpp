#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Matrice.hpp"
#include "Array2D.hpp"
#include "Lapack.hpp"
#include "Rand.hpp"


/*!Class that contains the information on the state
 * 
 * 
*/
template<typename T>
class System{
	public:
		/*!*/
		System();
		~System();

		unsigned int N_spin, N_m, N_site;
		
		void print();
		void init(unsigned int N_spin_, unsigned int N_m_, Matrice<double> const& H_, Array2D<unsigned int> const& sts_, Matrice<T> const& EVec, unsigned int thread);

		void update();
		void swap();

		void swap(unsigned int const& a, unsigned int const& b);
		//{Description
		/*!Computes the ratio of the two determinants related to the current and next
		 * configuration
		 *
		 * - when one matrix is modified, two of its columns are exchanged and
		 * therefore a minus sign arises 
		 * - when two matrices are modified, one computes the ratio using a specialized
		 * formula 
		 */
		//}
		T ratio();

		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		double compute_energy();

	private:
		/*!Forbids copy constructor*/
		System(System const& S);
		/*!Forbids assigment operator*/
		System& operator=(System const& S);

		Rand* rnd;			//!< generator of random numbers 
		Matrice<T> *A;      //!< det(A) <=> <GS|a>
		Matrice<T> *Ainv;   //!< inverse of A
		Matrice<T> tmp_mat; //!< temporary matrix used to update Ainv
		T w[2];             //!< determinant ratios : <GS|a>/<GS|b>
		unsigned int *s;    //!< s[i] = j : ??? ???
		unsigned int *wis;  //!< wis[i] = j : on ith site there is the j particle
		unsigned int mc[2]; //!< matrices (colors) that are modified 
		unsigned int cc[2]; //!< column's matrices (~band) that are exchanged 
		Matrice<double> H;	//!< Hamiltonian
		Array2D<unsigned int> sts; //!< sts(i,0) is a neighboor of sts(i,1)
};

/*constructors and destructor*/
/*{*/
template<typename T>
System<T>::System():
	N_spin(0),
	N_m(0),
	N_site(0),
	rnd(NULL),
	A(NULL),
	Ainv(NULL),
	tmp_mat(0,0.0),
	s(NULL),
	wis(NULL),
	H(0,0),
	sts(0,0)
{ }

template<typename T>
System<T>::~System(){
	delete[] s;
	delete[] wis;
	delete[] A;
	delete[] Ainv;
	delete rnd;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename T> 
T System<T>::ratio(){
	if(mc[0] == mc[1]){
		return -1;
	} else {
		w[0] = 0.0;
		w[1] = 0.0;
		for(unsigned int i(0);i<N_m;i++){
			w[0] += Ainv[mc[0]](cc[0],i)*A[mc[1]](i,cc[1]);
			w[1] += Ainv[mc[1]](cc[1],i)*A[mc[0]](i,cc[0]);
		}
		return w[0]*w[1]; 
	}
}

template<typename T>
double System<T>::compute_energy(){
	double E_step(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		E_step += real(ratio() * H(sts(j,0),sts(j,1)));
	}
	return E_step;
}

template<typename T>
void System<T>::print(){
	std::cout<<"{";
	for(unsigned int i(0); i<N_spin; i++){
		std::cout<<"{ ";
		for(unsigned int j(0); j<N_m; j++){
			std::cout<<s[i*N_m+j]<<",";
		}
		std::cout<<"},";
	}
	std::cout<<"}"<<std::endl;
	//for(unsigned int i(0); i<N_site; i++){
	//std::cout<<wis[i]<<" ";
	//}
	//std::cout<<std::endl;
	for(unsigned int i(0);i<2;i++){
		std::cout<<mc[i]<<" "<< cc[i]<<" "<<w[i]<<std::endl;
	}
	for(unsigned int i(0);i<N_spin;i++){
		std::cout<<std::endl;
		std::cout<<Ainv[i] * A[i]<<std::endl;
	}
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
	tmp_mat = Matrice<T>(N_m);
	A = new Matrice<T>[N_spin];
	Ainv = new Matrice<T>[N_spin];
	s = new unsigned int[N_site];
	wis = new unsigned int[N_site];
	rnd = new Rand(10,thread);

	unsigned int site(0);
	unsigned int N_as(N_site);
	unsigned int* available_sites(new unsigned int[N_site]);

	for(unsigned int i(0); i < N_as; i++){
		available_sites[i]  = i;	
	}
 
	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrice<T> (N_m);
		Ainv[i] = Matrice<T> (N_m);
		for(unsigned int j(0); j < N_m; j++){
			site = rnd->get(N_as);
			s[j+i*N_m] = available_sites[site];
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
	}

	for(unsigned int i(0);i<2;i++){
		w[i] = 0;
		cc[i] = 0;
		mc[i] = 0;
	}

	delete[] available_sites;
}

template<typename T>
void System<T>::swap() {
	mc[0] = rnd->get(N_spin);
	mc[1] = rnd->get(N_spin);
	while(mc[0]==mc[1]){
		mc[1] = rnd->get(N_spin);
	}
	cc[0] = rnd->get(N_m);
	cc[1] = rnd->get(N_m);
}

template<typename T>
void System<T>::swap(unsigned int const& a, unsigned int const& b) {
	mc[0] = wis[a] / N_m; //gives the color of particle on site a
	cc[0] = wis[a] % N_m; //gives the band of particle on site a
	mc[1] = wis[b] / N_m; //gives the color of particle on site b
	cc[1] = wis[b] % N_m; //gives the band of particle on site b
}

template<typename T>
void System<T>::update(){
	unsigned int a(0),b(0);
	a = mc[0]*N_m+cc[0];
	b = mc[1]*N_m+cc[1];

	unsigned int s_tmp(0);
	s_tmp = s[b];
	s[b] = s[a];
	s[a] = s_tmp;

	a = s[a];
	b = s[b];

	s_tmp = wis[b];
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
	T tmp_start(0.0);
	for(unsigned int m(0);m<2;m++){
		for(unsigned int i(0);i<N_m;i++){
			if(cc[m] == i){ tmp_start = -1.0; }
			else { tmp_start = 0.0; }
			for(unsigned int j(0);j<N_m;j++){
				tmp = tmp_start;
				for(unsigned int k(0);k<N_m;k++){
					tmp += Ainv[mc[m]](i,k)*A[mc[m]](k,cc[m]);
				}
				tmp_mat(i,j) = tmp*Ainv[mc[m]](cc[m],j)/w[m];
			}
		}
		Ainv[mc[m]] -= tmp_mat;
	}
}
/*}*/

/*double real(T)*/
/*{*/
template<typename T>
double real(T x);

template<>
inline double real(double x){
	return x;
}

template<>
inline double real(std::complex<double> x){
	return std::real(x);
}
/*}*/

#endif
