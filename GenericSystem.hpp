#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "PSTricks.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"
#include "IOSystem.hpp"
#include "Combination.hpp"

/*!Abstract class that can produce any kind of system*/
template<typename Type>
class GenericSystem:public Bosonic<Type>, public Fermionic<Type>, public IOSystem{
	public:
		/*{Description*/
		/*!Constructor requiring the coordination number and the name of the
		 * system, all other parameters of System have already been set by the
		 * most derived class (multiple inheritance)*/
		/*}*/
		GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename); 
		/*!Destructor*/
		virtual ~GenericSystem(){}

		/*{Description*/
		/*!Returns true if the system if bosonic 
		 * \warning not correctly implemented*/
		/*}*/
		bool is_bosonic() const { return false; }
		/*{Description*/
		/*!Saves ref_, N_, m_, n_, M_ and bc_ in jd_write_. As the method is
		 * virtual, a call on this method will call first child::save() const
		 * if it exists*/
		/*}*/
		virtual void save() const;
		virtual void create() = 0;
		virtual void check() = 0;

	protected:
		unsigned int const spuc_;//!< site per unit cell
		unsigned int const z_;	//!< coordination number
		RST system_info_;		//!< store information about the system

		/*{Description*/
		/*!Returns the neighbours of site i. 
		 *
		 * This pure virtual method must be defined here because it is needed
		 * by GenericSystem<Type>::compute_links()
		 */
		/*}*/
		virtual Matrix<int> get_neighbourg(unsigned int i) const = 0;
		/*!Computes the array of links between neighbouring sites*/
		void compute_links();
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename): 
	IOSystem(filename),
	spuc_(spuc),
	z_(z)
{ 
	filename_ += "-N" + tostring(this->N_);
	path_ += "N" + tostring(this->N_);
	filename_ += "-m" + tostring(this->m_);
	path_ += "/m" + tostring(this->m_);
	filename_ += "-n" + tostring(this->n_);
	path_ += "/n" + tostring(this->n_);
	filename_ += "-M";
	path_ += "/M";
	for(unsigned int i(0);i<this->M_.size();i++){
		filename_  += "-" + tostring(this->M_(i));
		path_ +=  "-"+tostring(this->M_(i));
	}
	switch(this->bc_){
		case -1:{filename_ += "-A"; path_ += "/A/"; }break;
		case 0: {filename_ += "-O"; path_ += "/O/"; }break;
		case 1: {filename_ += "-P"; path_ += "/P/"; }break;
		default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
	}
	path_ += tostring(this->ref_(0))+tostring(this->ref_(1))+tostring(this->ref_(2))+"/";
}

template<typename Type>
void GenericSystem<Type>::compute_links(){
	unsigned int z_link(0);
	unsigned int incr(0);
	if(z_%2){
		if(this->n_%2){ 
			std::cerr<<"void GenericSystem<Type>::compute_links() : can't compute links_"<<std::endl;
		} else {
			z_link = z_;
			incr = 2;
		}
	} else {
		z_link = z_/2;
		incr = 1;
	}
	unsigned int k(0);
	Matrix<int> nb;
	for(unsigned int i(0);i<this->n_;i+=incr){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<z_link;j++){
			if(nb(j,1)!=0){ k++; }
		}
	}
	this->links_.set(k,2);
	k=0;
	for(unsigned int i(0);i<this->n_;i+=incr){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<z_link;j++){
			if(nb(j,1)!=0){
				this->links_(k,0) = i;
				this->links_(k,1) = nb(j,0);
				k++;
			}
		}
	}
}

template<typename Type>
void GenericSystem<Type>::save() const {
	jd_write_->add_header()->add(system_info_.get());
	jd_write_->write("ref (type of wavefunction)",this->ref_);
	jd_write_->write("N (N of SU(N))",this->N_);
	jd_write_->write("m (# of particles per site)",this->m_);
	jd_write_->write("n (# of site)",this->n_);
	jd_write_->write("M (# of particles of each color, "+tostring(this->M_(0))+")",this->M_);
	jd_write_->write("bc (boundary condition)",this->bc_);
}
#endif
