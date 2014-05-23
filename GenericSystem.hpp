#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "Lapack.hpp"
#include "PSTricks.hpp"
#include "System.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.  */
template<typename Type>
class GenericSystem:public System, public Bosonic<Type>, public Fermionic<Type>{
	public:
		/*!Constructor N, n, m and z is the coordination number*/
		GenericSystem(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref, unsigned int const& z, std::string const& filename); 
		/*Simple destructor*/
		virtual ~GenericSystem();

		Matrix<unsigned int> get_links() const { return links_;}
		unsigned int get_num_links() const { return links_.row();}
		virtual	std::string get_filename() const = 0;

		virtual void create(double const& param, unsigned int const& type) = 0;
		virtual void save(IOFiles& w) const;
		virtual void check() = 0;

		bool is_bosonic() { return false; }
		bool is_degenerate() const { return degenerate_; }

	protected:
		unsigned int const z_;		//!< coordination number
		bool degenerate_;			//!< no degeneracy at the fermi level
		std::string filename_;		//!< filename of the System
		RST rst_;

		/*!return the neighbours of site i*/
		virtual Matrix<int> get_neighbourg(unsigned int i) const = 0;
		/*!compute the array of pairs of swapping sites*/
		void compute_links();
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref, unsigned int const& z, std::string const& filename): 
	System(N,n,m,bc,ref),
	z_(z),
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
void GenericSystem<Type>::compute_links(){
	unsigned int k(0);
	Matrix<int> nb;
	for(unsigned int i(0); i<n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			if(nb(j,1)!=0){ k++; }
		}
	}
	this->links_.set(k,2);
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
	Lapack<Type> ES(&this->T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	if(std::abs(EVal(M_) - EVal(M_-1))<1e-12){
		std::cerr<<"Fermi level degenerate"<<std::endl;
		degenerate_ = true;
	} else {
		degenerate_ = false;
	}
}

template<typename Type>
void GenericSystem<Type>::save(IOFiles& w) const {
	w.add_to_header(rst_.get());
	w("ref (type of wavefunction)",ref_);
	w("N (N of SU(N))",N_);
	w("m (# of particles per site)",m_);
	w("n (# of site)",n_);
	w("bc (boundary condition)",bc_);
}
#endif
