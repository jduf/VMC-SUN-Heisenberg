#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include <set>
#include "PSTricks.hpp"
#include "Gnuplot.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"
#include "IOSystem.hpp"
#include "Rand.hpp"
#include "Fit.hpp"

/*{*//*!Abstract class used by CreateSystem to produce any wavefunction
	   (Fermionic or Bosonic)

	   This class makes the connection between CreateSystem and more
	   specialized System like ChainPolymerized, SquareChiral... This is the
	   reason why it is an abstract class with many (pure) virtual methods.

	   Any instance of a child of this class will, at least, contain everything
	   that is required to compute the variational energy of the child related
	   trial wavefunction.*//*}*/
template<typename Type>
class GenericSystem: public Bosonic<Type>, public Fermionic<Type>, public IOSystem{
	public:
		/*{*//*!Constructor requiring only local parameters.
			   All other parameters of System have already been set by the most
			   derived class (multiple inheritance)*//*}*/
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

		/*!Sets bond energy(>=1), long range correlations(>=2)*/
		virtual void create_obs(unsigned int const& which_obs) = 0;

	protected:
		Matrix<Type> H_;					//!< matrix used to get the band structure
		Matrix<std::complex<double> > evec_;//!< eigenvectors of H+Tx+Ty
		unsigned int const spuc_;			//!< site per unit cell
		unsigned int const z_;	 			//!< coordination number
		double const eq_prec_ = 1e-12;		//!< precision for equality (important for matching position in lattice)

		/*{*//*!Returns the neighbours of site i.
			   This pure virtual method must be defined here because it is
			   needed by GenericSystem<Type>::set_nn_links() *//*}*/
		virtual Matrix<int> get_neighbourg(unsigned int const& i) const = 0;
		/*{*//*! Computes the array of links between neighbouring sites.
			   The argument l gives the number of links that need to be
			   computed for the site i%l.size() *//*}*/
		void set_nn_links(Vector<unsigned int> const& l);
		/*!Returns the index of the site i in the unit cell*/
		virtual unsigned int get_site_in_unit_cell(unsigned int const& i) const = 0;

		/*!Diagonalize the trial Hamiltonian H_*/
		void diagonalize(bool simple);
		/*!Evaluates the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);

	private:
		std::vector<std::string> generate_names() const;

		/*!Diagonalizes H_*/
		bool simple_diagonalization();
		/*!Diagonalizes H_+Translation operators => compute the band structure E(p)*/
		virtual bool full_diagonalization() = 0;
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	IOSystem(filename,generate_names()),
	spuc_(spuc),
	z_(z)
{
	if(this->status_ <= 4){
		if(this->status_ < 4){ std::cout<<__PRETTY_FUNCTION__<<" : should never happen, if it never does, might get rid of this extra test "<<this->status_<<std::endl;  }
		if(this->n_%this->spuc_){ std::cerr<<__PRETTY_FUNCTION__<<" : the number of sites ("<<this->n_<<") is not commensurate with the unit cell ("<<this->spuc_<<")"<<std::endl; }
		else { this->status_ = 3; }
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
		unsigned int l_tmp;
		Matrix<int> nb;
		if(!this->bc_){
			for(unsigned int i(0);i<this->n_;i++){
				l_tmp =l(i%l.size());
				if(l_tmp){
					nb = get_neighbourg(i);
					for(unsigned int j(0);j<l_tmp;j++){
						if(this->bc_ || nb(j,1)==0){ k++; }
					}
				}
			}
		} else { k = this->n_*this->z_/2; }
		Matrix<int> tmp(k,7);
		k=0;
		typedef std::tuple<unsigned int,unsigned int,unsigned int> ui_tuple ;
		std::set<ui_tuple> unit_cell_links;
		for(unsigned int i(0);i<this->n_;i++){
			l_tmp =l(i%l.size());
			if(l_tmp){
				nb = get_neighbourg(i);
				for(unsigned int j(0);j<l_tmp;j++){
					if(this->bc_ || nb(j,2)==0 ){
						tmp(k,0) = i;		//! site i
						tmp(k,1) = nb(j,0); //! site j
						tmp(k,3) = nb(j,1);	//! direction of vector linking i->j
						tmp(k,4) = nb(j,2); //! boundary condition test
						tmp(k,5) = get_site_in_unit_cell(i);
						tmp(k,6) = get_site_in_unit_cell(nb(j,0));
						unit_cell_links.insert(ui_tuple(tmp(k,3),tmp(k,5),tmp(k,6)));
						k++;
					}
				}
			}
		}
		if(unit_cell_links.size() != z_*spuc_/2){ std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number of links in the unit cell : "<<unit_cell_links.size() <<std::endl;; }
		else {
			/*!This value has nothing to do with the index of the bond l having
			 * a coupling J_(l). This value is only the index of a given bond
			 * in the unit cell. It means that J_(tmp(k,2)) would be completely
			 * absurd and wrong. The l-th link involved in the computation of
			 * E, has a bond energy of J_(l) and corresponds to the tmp(l,2)-th
			 * link the unit cell*/
			for(k=0;k<tmp.row();k++){
				tmp(k,2) = std::distance(unit_cell_links.begin(),unit_cell_links.find(ui_tuple(tmp(k,3),tmp(k,5),tmp(k,6))));
			}
			this->obs_.push_back(Observable("Energy per site",0,1,tmp));
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number of link"<<std::endl; }
	if(!this->bc_){ std::cerr<<__PRETTY_FUNCTION__<<" : open boundary condition could be problematic when nb(j,1)=0 and l(j) != 0"<<std::endl; }
}

template<typename Type>
std::vector<std::string> GenericSystem<Type>::generate_names() const {
	std::vector<std::string> parameter_names;
	parameter_names.push_back("N" + my::tostring(this->N_));
	parameter_names.push_back("m" + my::tostring(this->m_));
	parameter_names.push_back("n" + my::tostring(this->n_));
	std::string tmp("M");
	for(unsigned int i(0);i<this->M_.size();i++){ tmp  += "_" + my::tostring(this->M_(i)); }
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

template<typename Type>
void GenericSystem<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_ = 1; } }
	else { if(full_diagonalization()){ this->status_ = 1; } }
}

template<typename Type>
bool GenericSystem<Type>::simple_diagonalization(){
	Vector<double> eval;
	Lapack<Type>(H_,false,(this->ref_(1)==1?'S':'H')).eigensystem(eval,true);
	for(unsigned int c(0);c<this->N_;c++){
		if(my::are_equal(eval(this->M_(c)),eval(this->M_(c)-1),1e-12)){
			std::cerr<<__PRETTY_FUNCTION__<<" : degenerate at the Fermi level"<<std::endl;
			return false;
		}
	}
	return true;
}

template<typename Type>
std::complex<double> GenericSystem<Type>::projection(Matrix<Type> const& O, unsigned int const& idx){
	std::complex<double> tmp;
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		tmp = 0.0;
		for(unsigned int j(0);j<O.col();j++){ tmp += O(i,j)*evec_(j,idx); }
		out += std::conj(evec_(i,idx))*tmp;
	}
	return out;
}
#endif
