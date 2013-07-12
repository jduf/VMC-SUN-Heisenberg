#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Parseur.hpp"
#include "Write.hpp"
#include "Read.hpp"
#include "Lapack.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
 *
 *  
 *
*/
template<typename Type>
class CreateSystem{
	public:
		CreateSystem(Parseur& P,unsigned int N_n);
		~CreateSystem();

	protected:
		unsigned int const N_m, N_spin, N_site;//!<
		double bc;//!<
		Matrix<unsigned int> sts;//!< 
		Matrix<double> H;//!< hopping matrix
		Matrix<Type> T;//!< eigenvectors matrix (transfer matrix)
		bool successful;//!<
		bool is_complex;//!<
		RST rst;//!< will be added before the values in the header
		Write* w;//!<output file that contains all the informations of the system
		std::string filename; //!<
		char mat_type;//!<

		/*!Compute the eigenvectors from the mean field hamiltonian*/
		void compute_EVec();
		/*!Compute the array of pairs of swapping sites*/
		void compute_sts();

		void save();
};

template<typename Type>
CreateSystem<Type>::CreateSystem(Parseur& P, unsigned int N_n):
	N_m(P.get<unsigned int>("N_m")),
	N_spin(P.get<unsigned int>("N_spin")), 
	N_site(N_spin*N_m),
	bc(0),
	sts(N_spin*N_m*N_n/2,2),
	H(N_spin*N_m,N_spin*N_m,0.0),
	T(N_spin*N_m,N_spin*N_m,0.0),
	successful(false),
	is_complex(false),
	rst(),
	w(NULL),
	filename(""),
	mat_type('U')
{ }

template<typename Type>
CreateSystem<Type>::~CreateSystem(){
	if(w){ delete w; }
}

template<typename Type>
void CreateSystem<Type>::compute_sts(){
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

template<typename Type>
void CreateSystem<Type>::compute_EVec(){
	successful = false;
	Lapack<Type> ES(&T,false, mat_type);
	Matrix<double> EVal;
	ES.eigensystem(EVal);
	if(std::abs(EVal(N_m) - EVal(N_m-1))>1e-10){ successful = true; }
}

template<typename Type>
void CreateSystem<Type>::save(){
	if(successful){
		w = new Write(filename+".jdbin");
		rst.title("Input values","~");
		w->set_header(rst.get());
		(*w)("is_complex",is_complex);
		(*w)("N_spin",N_spin);
		(*w)("N_m",N_m);
		(*w)("sts",sts);
		(*w)("T",T);
		(*w)("bc",bc);
	}
}
#endif
