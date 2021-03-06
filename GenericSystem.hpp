#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include <set>
#include "PSTricks.hpp"
#include "Gnuplot.hpp"
#include "Bosonic.hpp"
#include "Fermionic.hpp"
#include "IOSystem.hpp"
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
class GenericSystem:public Bosonic<Type>, public Fermionic<Type>, public IOSystem{
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

		/*!Creates observable (bond energy, long range correlations,...)*/
		void create_obs(unsigned int const& which_obs);
		/*!Substitute method for wavefunctions without parameters to save*/
		virtual void save_param(IOFiles& w) const;
		virtual void create() = 0;
		virtual void check() = 0;
		virtual void display_results() = 0;

	protected:
		std::string const wf_name_ = "";	//!< the name of the wave function
		Matrix<std::complex<double> > evec_;//!< eigenvectors of H+Tx+Ty
		Matrix<Type> H_;					//!< matrix used to get the band structure
		unsigned int const spuc_;			//!< site per unit cell
		unsigned int const z_;	 			//!< coordination number
		double const eq_prec_ = 1e-12;		//!< precision for equality (important for matching position in lattice)

		/*{*//*!Creates the energy observable.
			   Computes the array of links between neighbouring sites. The
			   argument l gives the number of links that need to be computed
			   for the site i%l.size(). Once this array is computed, it is set
			   to the energy observable *//*}*/
		void energy_obs(Vector<unsigned int> const& l);
		/*!Create the bond energy observables*/
		virtual void bond_energy_obs();
		/*!Create the long range correlation observables*/
		virtual void long_range_correlations_obs() = 0;
		/*!Create the color occupation observables*/
		virtual void color_occupation_obs();
		/*!Create the observable H*H allowing the computation of the energy variance*/
		void energy_variance_obs();

		/*{*//*!Returns the neighbours of site i.
			   This pure virtual method must be defined here because it is
			   needed by GenericSystem<Type>::set_nn_links() *//*}*/
		virtual Matrix<int> get_neighbourg(unsigned int const& i) const = 0;
		/*!Returns the index of the site i in the unit cell*/
		virtual unsigned int site_index_to_unit_cell_index(unsigned int const& i) const = 0;

		/*!Diagonalize the trial Hamiltonian H_*/
		void diagonalize(bool simple);
		/*!Evaluates the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);

		/*!Used in VMCSystematic to plot E(param) with dim(param)=1*/
		virtual std::string extract_level_9();
		/*!Extrapolates the energy in the thermodynamical limit*/
		virtual std::string extract_level_2();
		/*!Gets the parameter for the fit of the thermodynamical limit*/
		virtual void param_fit_therm_limit(std::string& f, std::string& param, std::string& range);

		/*!Write the title, the command line and link the lattice in the rst file*/
		void rst_file_set_default_info(std::string const& param, std::string const& title, std::string const& replace = "");

	private:
		/*!Diagonalizes H_*/
		bool simple_diagonalization();
		/*!Diagonalizes H_+Translation operators (compute the band structure)*/
		virtual bool full_diagonalization() = 0;
};

template<typename Type>
GenericSystem<Type>::GenericSystem(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	IOSystem(filename,this->N_,this->m_,this->n_,this->M_,this->bc_,this->ref_),
	wf_name_(filename),
	spuc_(spuc),
	z_(z)
{
	if(this->status_ <= 4){
		if(this->n_%this->spuc_){ std::cerr<<__PRETTY_FUNCTION__<<" : n="<<this->n_<<" is not commensurate with the unit cell size="<<this->spuc_<<std::endl; }
		else { this->status_ = 3; }
	}
}

/*{public methods*/
template<typename Type>
void GenericSystem<Type>::create_obs(unsigned int const& which_obs){
	switch(which_obs){
		/*As the energy_variance_obs observable is certainly not correctly
		 * defined, one needs to explicitly ask for its creation, i.e. if
		 * which_obs==0, it will not create this observable because of the -1
		 * in the for loop*/
		case 0: { for(unsigned int i(1);i<Observable::number_of_observables_defined-1;i++){ create_obs(i); } }break;
		case 1: { bond_energy_obs(); }break;
		case 2: { long_range_correlations_obs(); }break;
		case 3: { color_occupation_obs(); }break;
		case 4: { energy_variance_obs(); }break;
		default:{
					std::cerr<<__PRETTY_FUNCTION__<<" : unknown observable "<<which_obs<<std::endl;
					std::cerr<<"Available observables are :"<<std::endl;
					std::cerr<<" + Bond energy            : 1"<<std::endl;
					std::cerr<<" + Long range correlations: 2"<<std::endl;
					std::cerr<<" + Color occupation       : 3"<<std::endl;
					std::cerr<<" + H*H (energy variance)  : 4"<<std::endl;
				}
	}
}

template<typename Type>
void GenericSystem<Type>::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>();
		w.add_to_header()->title("No parameter",'<');
		w.add_to_header()->add(system_info_.get());
	}
}
/*}*/

/*{protected methods*/
template<typename Type>
void GenericSystem<Type>::energy_obs(Vector<unsigned int> const& l){
	bool conflict(false);
	for(auto const& o:this->obs_){ if(o.get_type() == 0){ conflict = true; } }
	if(conflict){ std::cerr<<__PRETTY_FUNCTION__<<" : energy observable already defined"<<std::endl; }
	else if(2*l.sum()==l.size()*z_){
		unsigned int n_links(0);
		unsigned int l_tmp;
		Matrix<int> nb;
		if(!this->bc_){
			/*!Need to count the number of links in the cluster when using open
			 * boundary condition*/
			for(unsigned int i(0);i<this->n_;i++){
				l_tmp = l(i%l.size());
				if(l_tmp){
					nb = get_neighbourg(i);
					for(unsigned int j(0);j<l_tmp;j++){
						if(!nb(j,2)){ n_links++; }
					}
				}
			}
		} else { n_links = this->n_*z_/2; }
		Matrix<int> tmp(n_links,7);
		typedef std::tuple<unsigned int,unsigned int,unsigned int> ui_tuple;
		std::set<ui_tuple> unit_cell_links;
		unsigned int k(0);
		for(unsigned int i(0);i<this->n_;i++){
			l_tmp =l(i%l.size());
			if(l_tmp){
				nb = get_neighbourg(i);
				for(unsigned int j(0);j<l_tmp;j++){
					if(this->bc_ || !nb(j,2)){
						tmp(k,0) = i;		//! site i
						tmp(k,1) = nb(j,0); //! site j
						tmp(k,3) = nb(j,1);	//! direction of vector linking i->j
						tmp(k,4) = nb(j,2); //! boundary condition test
						tmp(k,5) = site_index_to_unit_cell_index(i);
						tmp(k,6) = site_index_to_unit_cell_index(nb(j,0));
						unit_cell_links.insert(ui_tuple(tmp(k,3),tmp(k,5),tmp(k,6)));
						k++;
					}
				}
			}
		}
		/*!The second condition is to avoid a problem with open boundary
		 * condition and cluster of size n_=spuc_*/
		if(unit_cell_links.size() != z_*spuc_/2 && spuc_ != this->n_){
			std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number of links in the unit cell ("<<unit_cell_links.size()<<" insead of "<<z_*spuc_/2<<"),  they are :"<<std::endl;
			for(auto const& uic:unit_cell_links){
				std::cout<<std::get<0>(uic)<<" "<<std::get<1>(uic)<<" "<<std::get<2>(uic)<<std::endl;
			}
		} else {
			/*!This value has nothing to do with the index of the bond l having
			 * a coupling J_(l). This value is only the index of a given bond
			 * in the unit cell. It means that J_(tmp(i,2)) would be completely
			 * absurd and wrong. The l-th link involved in the computation of
			 * E, has a bond energy of J_(l) and corresponds to the tmp(l,2)-th
			 * link the unit cell*/
			for(unsigned int i(0);i<tmp.row();i++){
				tmp(i,2) = std::distance(unit_cell_links.begin(),unit_cell_links.find(ui_tuple(tmp(i,3),tmp(i,5),tmp(i,6))));
			}
			this->obs_.push_back(Observable("Energy per site",0,1,tmp));
			this->ref_(4) = 0;
			this->status_ = 2;
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number of link"<<std::endl; }
	if(!this->bc_){ std::cerr<<__PRETTY_FUNCTION__<<" : open boundary condition could be problematic when nb(j,1)=0 and l(j) != 0"<<std::endl; }
}

template<typename Type>
void GenericSystem<Type>::bond_energy_obs(){
	this->obs_.push_back(Observable("Bond energy",1,z_*this->spuc_/2,this->obs_[0].nlinks()));
	this->obs_.back().remove_links();
}

template<typename Type>
void GenericSystem<Type>::color_occupation_obs(){
	Matrix<int> tmp(this->n_,this->N_);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->N_;j++){
			tmp(i,j) = site_index_to_unit_cell_index(i)*this->N_+j;
		}
	}
	this->obs_.push_back(Observable("Color occupation",3,this->N_*this->spuc_,tmp,this->n_/this->spuc_));
}

template<typename Type>
void GenericSystem<Type>::energy_variance_obs(){
	this->obs_.push_back(Observable("H*H (energy variance)",4,1,0));
}

template<typename Type>
void GenericSystem<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_ = 1; } }
	else { if(full_diagonalization()){ this->status_ = 1; } }
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

template<typename Type>
std::string GenericSystem<Type>::extract_level_9(){
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_);
	gp+="plot '"+this->filename_+".dat' u 1:2:3 notitle";
	gp.save_file();
	gp.create_image(true,"png");
	return this->filename_;
}

template<typename Type>
std::string GenericSystem<Type>::extract_level_2(){
	std::string f("");
	std::string p("");
	std::string r("");
	param_fit_therm_limit(f,p,r);
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_);
	gp+=f;
	gp+="set fit quiet";
	gp+="fit " + r + " f(x) '"+this->filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) yerror via "+p;
	gp+="set print \'../"+this->sim_.substr(0,this->sim_.size()-1)+".dat\' append";
	gp+="print " + my::tostring(this->N_) + "," + my::tostring(this->m_) + "," + p[0];
	gp.range("x","0","");
	gp.label("x","$\\frac{ 1}{n}$");
	gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
	gp+="plot '"+this->filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) w e notitle,\\";
	gp+="     " + r + " f(x) t sprintf('%f'," + p[0] + ")";
	gp.save_file();
	gp.create_image(true,"png");

	return this->filename_;
}

template<typename Type>
void GenericSystem<Type>::param_fit_therm_limit(std::string& f, std::string& param, std::string& range){
	f = "f(x) = a+b*x";
	param = "a,b";
	range = "";
}

template<typename Type>
void GenericSystem<Type>::rst_file_set_default_info(std::string const& param, std::string const& title, std::string const& replace){
	if(this->rst_file_){
		std::string cmd_name("./mc");
		cmd_name+= " -s:wf "+this->wf_name_+ " " + param;
		cmd_name+= " -u:N " + my::tostring(this->N_);
		cmd_name+= " -u:m " + my::tostring(this->m_);
		cmd_name+= " -u:n " + my::tostring(this->n_);
		cmd_name+= " -i:bc "+ my::tostring(this->bc_);
		cmd_name+= " -d -u:tmax 10";

		if(!replace_title_with_link_in_rst_){
			this->rst_file_->title(title,'-');
		} else {
			this->rst_file_->title("|"+replace+"|_",'-');
			this->rst_file_->replace(replace,title);
		}
		this->rst_file_->change_text_onclick("run command",cmd_name);
		this->rst_file_->figure(
				this->dir_+this->filename_+".png",
				RST::math("E="+my::tostring(this->obs_[0][0].get_x())+"\\pm"+my::tostring(this->obs_[0][0].get_dx())),
				RST::target(this->dir_+this->filename_+".pdf")+RST::width("800")
				);
	}
}
/*}*/

/*{private methods*/
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
/*}*/
#endif
