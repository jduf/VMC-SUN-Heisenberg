#include "Write.hpp"
#include "Array2D.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"

#include <complex>
#include <sstream>

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
 *
 *  
 *
*/
template<typename M>
class CreateState{
	public:
		CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename);
		~CreateState();

	private:
		unsigned int const N_m, N_n, N_spin, N_site;//!<
		Array2D<unsigned int> sts;//!< 
		Matrice<double> H;//!< hopping matrix
		Matrice<M> T;//!< eigenvectors matrix (transfer matrix)
		std::string filename; //!<
		bool is_complex;//!<
		char mat_type;//!<

		/*!Compute the hopping matrix*/
		void compute_H();
		/*!Compute the eigenvectors from the mean field hamiltonian*/
		void compute_EVec();
		/*!For every site, compute all its other sites that  */
		void compute_sts();
};

template<typename M>
CreateState<M>::CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename):
	N_m(N_m),
	N_n(N_n),
	N_spin(N_spin), 
	N_site(N_m*N_spin),
	sts(N_spin*N_m*N_n/2,2),
	H(N_spin*N_m),
	T(N_spin*N_m),
	filename(filename),
	is_complex(false),
	mat_type('U')
{
	compute_H();
	compute_sts();
	compute_EVec();
}


template<typename M>
CreateState<M>::~CreateState(){
	Write w(filename+".jdbin");
	w("complex",is_complex);
	w("N_spin",N_spin);
	w("N_m",N_m);
	w("N_n",N_n);
	w("sts",sts);
	w("H",H);
	w("T",T);
}

template<typename M>
void CreateState<M>::compute_sts(){
	unsigned int k(0);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( std::abs(H(i,j)) > 1e-4){
				sts(k,0) = i;
				sts(k,1) = j;
				k++;
			}
		}
	}
}

template<typename M>
void CreateState<M>::compute_EVec(){
	Lapack<M> ES(T.ptr(),N_site, mat_type);
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	if(std::abs(EVal(N_m) - EVal(N_m-1))<1e-10){
		filename += "-VPDNF";
	}
}
