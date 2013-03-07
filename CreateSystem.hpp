#include "Read.hpp"
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
		std::string filename; //!<
		std::string filename_comp;//!<
		bool is_complex;//!<
		char mat_type;//!<
		Write w;//!<
		unsigned int const N_m, N_spin, N_n, N_site;//!<
		Matrice<M> H;//!< hopping matrix
		Matrice<M> EVec;//!< eigenvectors matrix
		Array2D<unsigned int> sts;//!<

		/*!Compute the hopping matrix*/
		void compute_H();
		/*!Compute the eigenvectors from H*/
		void compute_EVec();
		/*!For every site, compute all its other sites that  */
		void compute_sts();
};

template<typename M>
CreateState<M>::CreateState(unsigned int N_m, unsigned int N_spin, unsigned int N_n, std::string filename):
	filename(filename),
	filename_comp(""),
	is_complex(false),
	mat_type('U'),
	w(),
	N_m(N_m),
	N_spin(N_spin), 
	N_n(N_n),
	N_site(N_m*N_spin),
	H(N_spin*N_m),
	EVec(N_spin*N_m),
	sts(N_spin*N_m*N_n/2,2)
{
	compute_H();
	compute_sts();
	compute_EVec();
}


template<typename M>
CreateState<M>::~CreateState(){
	w.open(filename+filename_comp+".jdbin");
	w("complex",is_complex);
	w("N_spin",N_spin);
	w("N_m",N_m);
	w("N_n",N_n);
	w("sts",sts);
	w("H",H);
	w("T",EVec);
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
	EVec = H;
	Lapack<M> ES(EVec.ptr(),N_site, mat_type);
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	if(std::abs(EVal(N_m) - EVal(N_m-1))<1e-10){
		filename_comp = "-VPDNF";
	}
}
