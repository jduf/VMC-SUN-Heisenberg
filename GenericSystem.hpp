#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "PSTricks.hpp"
#include "Gnuplot.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"
#include "IOSystem.hpp"
#include "Rand.hpp"
#include "Fit.hpp"

/*!Abstract class that can produce any kind of system*/
template<typename Type>
class GenericSystem: public Bosonic<Type>, public Fermionic<Type>, public IOSystem{
	public:
		/*{Description*/
		/*!Constructor requiring the coordination number and the name of the
		 * system, all other parameters of System have already been set by the
		 * most derived class (multiple inheritance)*/
		/*}*/
		GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename);
		/*!Default destructor*/
		virtual ~GenericSystem() = default;
		/*{Forbidden*/
		GenericSystem() = delete;
		GenericSystem(GenericSystem<Type> const&) = delete;
		GenericSystem(GenericSystem<Type>&&) = delete;
		GenericSystem& operator=(GenericSystem<Type> const&) = delete;
		/*}*/

		virtual void save_param(IOFiles& w) const;
		virtual void create() = 0;
		virtual void check() = 0;
		virtual void display_results() = 0;
		virtual void get_wf_symmetries(std::vector<Matrix<int> >& sym) const { (void)(sym); }

		/*!Sets the binning for E_(>=0), corr_(>=1), lr_corr(>=2)*/
		virtual void set_observables(int nobs) = 0;

	protected:
		unsigned int const spuc_;//!< site per unit cell
		unsigned int const z_;	 //!< coordination number
		RST system_info_;		 //!< store information about the system

		/*{Description*/
		/*!Returns the neighbours of site i.
		 *
		 * This pure virtual method must be defined here because it is needed
		 * by GenericSystem<Type>::set_nn_links()
		 */
		/*}*/
		virtual Matrix<int> get_neighbourg(unsigned int const& i) const = 0;
		/*{Description*/
		/*!Computes the array of links between neighbouring sites
		 * The argument l gives the number of links that need to be computed
		 * for the site i%l.size()*/
		/*}*/
		void set_nn_links(Vector<unsigned int> const& l);
		void check_lattice();

	private:
		std::vector<std::string> generate_names() const;
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	IOSystem(filename,generate_names()),
	spuc_(spuc),
	z_(z)
{
	if(this->n_%this->spuc_){
		this->status_= 4;
		std::cerr<<__PRETTY_FUNCTION__<<" : the number of sites ("<<this->n_<<") is not commensurate with the unit cell ("<<this->spuc_<<")"<<std::endl;
	} else {
		/*!Need to redefine the value of status to be 3 because in MCSim::save,
		 * CreateSystem uses a MCSystem in the constructor. As
		 * MCSystem::stats_=0 for successful simulation, the created System
		 * will also have System::status_=0 which will be problematic.*/
		this->status_=3;
	}
}

template<typename Type>
void GenericSystem<Type>::save_param(IOFiles& w) const {
	w.add_header()->add(system_info_.get());
}

template<typename Type>
void GenericSystem<Type>::set_nn_links(Vector<unsigned int> const& l){
	if(2*l.sum()==l.size()*z_){
		unsigned int k(0);
		Matrix<int> nb;
		for(unsigned int i(0);i<this->n_;i++){
			nb = get_neighbourg(i);
			for(unsigned int j(0);j<l(i%l.size());j++){
				if(nb(j,1)!=0){ k++; }
			}
		}
		Matrix<int> tmp(k,5);
		k=0;
		for(unsigned int i(0);i<this->n_;i++){
			nb = get_neighbourg(i);
			for(unsigned int j(0);j<l(i%l.size());j++){
				if(nb(j,1)!=0){
					tmp(k,0) = i;
					tmp(k,1) = nb(j,0);
					/*!set tmp(k,2) to -1 so that when one wants to measure the
					 * bond energy, one has to redefine the correct mapping
					 * between all bonds and the representative ones*/
					tmp(k,2) = -1;
					tmp(k,3) = j;
					tmp(k,4) = nb(j,1);
					k++;
				}
			}
		}
		this->obs_.push_back(Observable(tmp));
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number of link"; }
	if(!this->bc_){ std::cerr<<__PRETTY_FUNCTION__<<" : open boundary condition could be problematic when nb(j,1)=0 and l(j) != 0"<<std::endl; }
}

template<typename Type>
void GenericSystem<Type>::check_lattice(){
	Vector<int> dir;
	Rand<unsigned int> rnd(0,z_-1);
	unsigned int j;
	unsigned int d;
	unsigned int p0,p1;
	std::cout<<"check"<<std::endl;
	for(unsigned int iter(10);iter<1e4;iter*=5){
		dir.set(z_,0);
		p0=0;
		for(unsigned int i(0);i<iter;i++){
			d=rnd.get();
			p0=get_neighbourg(p0)(d,0);
			dir(d)++;
		}
		d=p1=j=0;
		while(j<iter){
			if(dir(d)){
				j++;
				dir(d)--;
				p1=get_neighbourg(p1)(d,0);
			} else { d++; }
		}
		if(p0 != p1){
			this->status_++;
			std::cerr<<__PRETTY_FUNCTION__<<" : no consistent enumeration"<<std::endl;
		}
	}
}

template<typename Type>
std::vector<std::string> GenericSystem<Type>::generate_names() const {
	std::vector<std::string> parameter_names;
	parameter_names.push_back("N" + my::tostring(this->N_));
	parameter_names.push_back("m" + my::tostring(this->m_));
	parameter_names.push_back("n" + my::tostring(this->n_));
	std::string tmp("M");
	for(unsigned int i(0);i<this->M_.size();i++){
		tmp  += "_" + my::tostring(this->M_(i));
	}
	parameter_names.push_back(tmp);
	switch(this->bc_){
		case -1:{ parameter_names.push_back("A"); }break;
		case 0: { parameter_names.push_back("O"); }break;
		case 1: { parameter_names.push_back("P"); }break;
	}
	parameter_names.push_back("Juniform");
	parameter_names.push_back(my::tostring(this->ref_(0))+my::tostring(this->ref_(1))+my::tostring(this->ref_(2)));
	return parameter_names;
}
#endif
