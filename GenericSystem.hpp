#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "PSTricks.hpp"
#include "Gnuplot.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"
#include "IOSystem.hpp"
#include "Rand.hpp"

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
		virtual void lattice(std::string const& path) = 0;
		virtual void get_wf_symmetries(std::vector<Matrix<int> >& sym) const { (void)(sym); }

	protected:
		unsigned int const spuc_;//!< site per unit cell
		unsigned int const z_;	 //!< coordination number
		RST system_info_;		 //!< store information about the system

		/*{Description*/
		/*!Returns the neighbours of site i. 
		 *
		 * This pure virtual method must be defined here because it is needed
		 * by GenericSystem<Type>::compute_links()
		 */
		/*}*/
		virtual Matrix<int> get_neighbourg(unsigned int const& i) const = 0;
		/*{Description*/
		/*!Computes the array of links between neighbouring sites
		 * The argument l gives the number of links that need to be computed
		 * for the site i%l.size()*/
		/*}*/
		void compute_links(Vector<unsigned int> const& l);
		void check_lattice();
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename): 
	IOSystem(filename,this->names()),
	spuc_(spuc),
	z_(z)
{
	if(this->n_%this->spuc_){
		this->status_++;
		std::cerr<<"GenericSystem<Type>::GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename) : the number of sites is not comensurate with the unit cell"<<std::endl;
	}
}

template<typename Type>
void GenericSystem<Type>::save_param(IOFiles& w) const {
	w.add_header()->add(system_info_.get());
}

template<typename Type>
void GenericSystem<Type>::compute_links(Vector<unsigned int> const& l){
	if(2*l.sum()==l.size()*z_){
		unsigned int k(0);
		Matrix<int> nb;
		for(unsigned int i(0);i<this->n_;i++){
			nb = get_neighbourg(i);
			for(unsigned int j(0);j<l(i%l.size());j++){
				if(nb(j,1)!=0){ k++; }
			}
		}
		this->links_.set(k,2);
		k=0;
		for(unsigned int i(0);i<this->n_;i++){
			nb = get_neighbourg(i);
			for(unsigned int j(0);j<l(i%l.size());j++){
				if(nb(j,1)!=0){
					this->links_(k,0) = i;
					this->links_(k,1) = nb(j,0);
					k++;
				}
			}
		}
	} else {
		std::cerr<<"void GenericSystem<Type>::compute_links(Vector<unsigned int> const& l) : incoherent number of link";
	}
	if(this->bc_==0){
		std::cerr<<"void GenericSystem<Type>::compute_links(Vector<unsigned int> const& l) : open boundary condition could be problematic when nb(j,1)=0 and l(j) != 0"<<std::endl;
	}
}

template<typename Type>
void GenericSystem<Type>::check_lattice(){
	Matrix<int> nb;
	Vector<int> dir;
	Rand<unsigned int> rnd(0,z_-1);
	unsigned int j;
	unsigned int d;
	unsigned int p0,p1;
	std::cout<<"check"<<std::endl;
	for(unsigned int iter(10);iter<1e4;iter*=5){
		dir.set(z_,0);
		d=0;
		p0=0;
		for(unsigned int i(0);i<iter;i++){
			d=rnd.get();
			nb=get_neighbourg(p0);
			p0 = nb(d,0);
			dir(d)++;
		}
		d=0;
		p1=0;
		j=0;
		std::cout<<dir<<std::endl;
		while(j<iter){
			if(dir(d) != 0){
				dir(d)--;
				j++;
				nb=get_neighbourg(p1);
				p1 = nb(d,0);
			} else { d++; }
		}
		if(p0 != p1){ 
			this->status_++; 
			std::cerr<<"void GenericSystem<Type>::check_lattice() const : no consistent enumeration"<<std::endl;
		}
	}
}
#endif
