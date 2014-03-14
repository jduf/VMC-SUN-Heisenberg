#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "Lapack.hpp"
#include "PSTricks.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
*/
template<typename Type>
class GenericSystem{
	public:
		/*!Constructor N, n, m and z is the coordination number*/
		GenericSystem(unsigned int N, unsigned int n, unsigned int m, int bc, unsigned int z, std::string filename); 
		/*Simple destructor*/
		virtual ~GenericSystem();

		Matrix<unsigned int> get_links() const { return links_;}
		Matrix<Type> get_EVec() const { return EVec_;}
		unsigned int get_num_links() const { return links_.row();}
		std::string get_filename() const { return filename_;}

		virtual unsigned int create(double param)=0;
		virtual void save(Write& w) const;
		virtual void check()=0;
		virtual void study(double E, double DeltaE, Vector<double> const& corr, std::string save_in) = 0;

	protected:
		unsigned int const n_;		//!< sites' number
		unsigned int const N_;		//!< colors' number
		unsigned int const m_;		//!< particles per site's number
		unsigned int const M_;		//!< particles' number of each color
		unsigned int const z_;		//!< coordination number
		int bc_;					//!< boundary condition
		Matrix<unsigned int> links_;//!< list of links
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		bool degenerate_;			//!< no degeneracy at the fermi level
		std::string filename_;		//!< filename of the System
		RST rst_;

		/*!return the neighbours of site i*/
		virtual Matrix<int> get_neighbourg(unsigned int i) const=0;
		/*!compute the array of pairs of swapping sites*/
		void compute_sts();
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int N, unsigned int n, unsigned int m, int bc, unsigned int z, std::string filename): 
	n_(n),
	N_(N), 
	m_(m),
	M_((m_*n_)/N_), 
	z_(z),
	bc_(bc),
	T_(n_,n_,0.0),
	EVec_(N_*n_,M_),
	degenerate_(false),
	filename_(filename)
{ 
	if(M_*N_ != m_*n_){ std::cerr<<"GenericSystem::GenericSystem(N,n,m,z,filename) : there is not an equal number of color"<<std::endl; }
	if(m_ > N_){ std::cerr<<"GenericSystem::GenericSystem(N,n,m,z,filename) : m>N is impossible"<<std::endl;} 
	filename_ += "-N" + tostring(N_);
	filename_ += "-m" + tostring(m_);
	filename_ += "-n" + tostring(n_);
	switch(bc_){
		case -1:{filename_ += "-A"; }break;
		case 0: {filename_ += "-O"; }break;
		case 1: {filename_ += "-P"; }break;
		default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
	}
}

template<typename Type>
GenericSystem<Type>::~GenericSystem(){}

template<typename Type>
void GenericSystem<Type>::compute_sts(){
	unsigned int k(0);
	Matrix<int> nb;
	for(unsigned int i(0); i<n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			if(nb(j,1)!=0){ k++; }
		}
	}
	links_.set(k,2);
	k=0;
	for(unsigned int i(0); i<n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			if(nb(j,1)!=0){
				links_(k,0) = i;
				links_(k,1) = nb(j,0);
				k++;
			}
		}
	}
}

template<typename Type>
void GenericSystem<Type>::diagonalize_T(char mat_type){
	Lapack<Type> ES(&T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	if(std::abs(EVal(M_) - EVal(M_-1))<1e-12){
		std::cerr<<"Fermi level degenerate"<<std::endl;
		degenerate_ = true;
	}
}

template<typename Type>
void GenericSystem<Type>::save(Write& w) const {
	w.add_to_header(rst_.get());
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("n (particles' number)",n_);
	w("bc (boundary condition)",bc_);
}
#endif
