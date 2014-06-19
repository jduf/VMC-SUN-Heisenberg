#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "PSTricks.hpp"
#include "BandStructure.hpp"
#include "RSTFile.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.  */
template<typename Type>
class GenericSystem:public Bosonic<Type>, public Fermionic<Type>{
	public:
		/*!Constructor N, n, m and z is the coordination number*/
		GenericSystem(unsigned int const& z, std::string const& filename); 
		/*Simple destructor*/
		virtual ~GenericSystem();

		unsigned int get_num_links() const { return this->links_.row();}
		virtual	std::string get_filename() const = 0;

		virtual void create(double const& param, unsigned int const& type) = 0;
		virtual void save(IOFiles& w) const;
		virtual void check() = 0;

		bool is_bosonic() { return false; }

	protected:
		unsigned int const z_;		//!< coordination number
		std::string filename_;		//!< filename of the System
		RST rst_;

		/*!return the neighbours of site i*/
		virtual Matrix<int> get_neighbourg(unsigned int i) const = 0;
		/*!compute the array of pairs of swapping sites*/
		void compute_links();
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& z, std::string const& filename): 
	z_(z),
	filename_(filename)
{ 
	std::cout<<"genericsystem"<<std::endl;
	if(this->M_*this->N_ != this->m_*this->n_){ std::cerr<<"GenericSystem::GenericSystem(N,n,m,z,filename) : there is not an equal number of color"<<std::endl; }
	if(this->m_>this->N_){ std::cerr<<"GenericSystem::GenericSystem(N,n,m,z,filename) : m>N is impossible"<<std::endl;} 
	filename_ += "-N" + tostring(this->N_);
	filename_ += "-m" + tostring(this->m_);
	filename_ += "-n" + tostring(this->n_);
	switch(this->bc_){
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
	for(unsigned int i(0); i<this->n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			if(nb(j,1)!=0){ k++; }
		}
	}
	this->links_.set(k,2);
	k=0;
	for(unsigned int i(0); i<this->n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			if(nb(j,1)!=0){
				this->links_(k,0) = i;
				this->links_(k,1) = nb(j,0);
				k++;
			}
		}
	}
}


template<typename Type>
void GenericSystem<Type>::save(IOFiles& w) const {
	w.add_to_header(rst_.get());
	w("ref (type of wavefunction)",this->ref_);
	w("N (N of SU(N))",this->N_);
	w("m (# of particles per site)",this->m_);
	w("n (# of site)",this->n_);
	w("bc (boundary condition)",this->bc_);
}
#endif
