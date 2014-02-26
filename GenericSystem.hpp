#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "Lapack.hpp"
#include "Gnuplot.hpp"
#include "PSTricks.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
*/
template<typename Type>
class GenericSystem{
	public:
		/*!Constructor N, n, m and z is the coordination number*/
		GenericSystem(unsigned int N, unsigned int n, unsigned int m, unsigned int z, std::string filename); 
		/*Simple destructor*/
		virtual ~GenericSystem();

		Matrix<unsigned int> get_sts() const { return sts_;}
		Matrix<Type> get_EVec() const { return EVec_;}
		unsigned int get_num_links() const { return n_*z_/2;}
		std::string get_filename() const { return filename_;}

		virtual void save(Write& w);
		virtual void create(double param)=0;

	protected:
		unsigned int const n_;		//!< sites' number
		unsigned int const N_;		//!< colors' number
		unsigned int const m_;		//!< particles per site's number
		unsigned int const M_;		//!< particles' number of each color
		unsigned int const z_;		//!< coordination number
		double bc_;					//!< boundary condition
		Matrix<unsigned int> sts_;	//!< list of connected sites
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		bool degenerate_;			//!< no degeneracy at the fermi level
		std::string filename_;		//!< filename of the System
		RST rst_;

		/*!return the neighbours of site i*/
		virtual Vector<unsigned int> get_neighbourg(unsigned int i)=0;
		/*!compute the array of pairs of swapping sites*/
		void compute_sts();
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int N, unsigned int n, unsigned int m, unsigned int z, std::string filename): 
	n_(n),
	N_(N), 
	m_(m),
	M_((m_*n_)/N_), 
	z_(z),
	bc_(1),
	sts_(n_*z_/2,2),
	T_(n_,n_,0.0),
	EVec_(N_*n_,M_),
	degenerate_(false),
	filename_(filename)
{ 
	if(M_*N_ != m_*n_){ std::cerr<<"GenericSystem::GenericSystem(N,n,m,z,filename) : there is not an equal number of color"<<std::endl; }
	if(m_ > N_){ std::cerr<<"GenericSystem::GenericSystem(N,n,m,z,filename) : m>N is impossible"<<std::endl;} 
	filename_ += "-N" + tostring(N_);
	filename_ += "-m" + tostring(m_);
	filename_ += "-S" + tostring(n_);
}

template<typename Type>
GenericSystem<Type>::~GenericSystem(){}

template<typename Type>
void GenericSystem<Type>::compute_sts(){
	unsigned int k(0);
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0); i<n_;i++){
		neighbourg = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			sts_(k,0) = i;
			sts_(k,1) = neighbourg(j);
			k++;
		}
	}
}

template<typename Type>
void GenericSystem<Type>::diagonalize_T(char mat_type){
	Lapack<Type> ES(&T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	if(std::abs(EVal(M_) - EVal(M_-1))<1e-12){ degenerate_ = true; }
}

template<typename Type>
void GenericSystem<Type>::save(Write& w){
	w.set_header(rst_.get());
	w("n (particles' number)",n_);
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("bc (boundary condition)",bc_);
}
#endif
