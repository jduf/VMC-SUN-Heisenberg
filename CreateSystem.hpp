#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Parseur.hpp"
#include "Lapack.hpp"
#include "Gnuplot.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
 *
 *  
 *
*/
template<typename Type>
class CreateSystem{
	public:
		/*!Parseur needs N and m, z is the coordination number*/
		CreateSystem(Parseur& P, unsigned int z, std::string filename); 
		/*Simple destructor*/
		~CreateSystem();

	protected:
		std::string wf_;			//!< type of wavefunction
		unsigned int const m_;		//!< number of unit cell
		unsigned int const N_;		//!< N of SU(N)
		unsigned int const n_;		//!< number of sites
		unsigned int const z_;		//!< coordination number
		double bc_;					//!< boundary condition
		Matrix<unsigned int> sts_;	//!< list of connected sites
		Matrix<int> H_;				//!< SU(N) Matrix
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		bool successful_;			//!< no degeneracy at the fermi level
		std::string filename_;

		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
		/*!compute the array of pairs of swapping sites*/
		void compute_sts();
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket);

		void save_band_structure(Vector<double> kx, Vector<double> ky, Vector<double> E);
};

template<typename Type>
CreateSystem<Type>::CreateSystem(Parseur& P, unsigned int z, std::string filename): 
	wf_(P.get<std::string>("wf")),
	m_(P.get<unsigned int>("m")),
	N_(P.get<unsigned int>("N")), 
	n_(N_*m_),
	z_(z),
	bc_(0),
	sts_(n_*z/2,2),
	H_(n_,n_,0.0),
	T_(n_,n_,0.0),
	EVec_(N_*n_,m_),
	successful_(false),
	filename_(filename)
{ }

template<typename Type>
CreateSystem<Type>::~CreateSystem(){ }

template<typename Type>
void CreateSystem<Type>::compute_sts(){
	unsigned int k(0);
	for(unsigned int i(0); i<n_;i++){
		for(unsigned int j(i+1); j<n_;j++){
			if ( H_(i,j) != 0){
				sts_(k,0) = i;
				sts_(k,1) = j;
				k++;
			}
		}
	}
}

template<typename Type>
void CreateSystem<Type>::diagonalize_T(char mat_type){
	Lapack<Type> ES(&T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	//std::cout<<EVal<<std::endl;
	if(std::abs(EVal(m_) - EVal(m_-1))>1e-10){ successful_ = true; }
}

template<typename Type>
std::complex<double> CreateSystem<Type>::projection(Matrix<Type> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket){
	Vector<std::complex<double> > tmp(O.row(),0.0);
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		for(unsigned int j(0);j<O.col();j++){
			tmp(i) += O(i,j)*base(j,ket);
		}
	}
	for(unsigned int i(0);i<O.row();i++){
		out += std::conj(base(i,bra))*tmp(i);
	}
	return out;
}

template<typename Type>
void CreateSystem<Type>::save_band_structure(Vector<double> kx, Vector<double> ky, Vector<double> E){
	Gnuplot gp(filename_+"-band-structure","splot");
	gp.save_data(filename_+"-spectrum",kx,ky,E);
	gp.add_plot_param(" ,\\\n");
	Vector<unsigned int> index(E.sort());
	gp.save_data(filename_+"-spectrum-sorted",kx.sort(index).range(0,m_),ky.sort(index).range(0,m_),E.range(0,m_));
}
#endif
