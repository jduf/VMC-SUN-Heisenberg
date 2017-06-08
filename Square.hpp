#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "System2D.hpp"

template<typename Type>
class Square: public System2D<Type>{
	public:
		/*!Constructor that organises the n=p^2+q^2 sites (p,q integer)*/
		Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Square()=0;

	protected:
		void init_lattice();
		/*Draw the lattice inside a PSTricks file*/
		void draw_lattice(bool const& only_unit_cell, bool const& silent, bool const& only_lattice, Vector<double> const& uc_shift, std::string param, std::string title);
		/*!To get the data for a plot of E(Jp) for all N,m,n,bc*/
		std::string extract_level_6();

	private:
		unsigned int p_;
		unsigned int q_;

		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& k, unsigned int const& spuc, unsigned int& ref3);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const;
};

/*{constructor*/
template<typename Type>
Square<Type>::Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry((this->ref_(4)?this->n_:0),this->N_/this->m_,spuc,this->ref_(3)),ab,spuc,4,4,filename)
{ this->filename_ += +"-p"+my::tostring(p_)+"-q"+my::tostring(q_); }

template<typename Type>
Square<Type>::~Square() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Square<Type>::init_lattice(){
	if(!this->ref_(4)){ this->status_ = 2; }
	else {
		if(this->dir_nn_){
			/*{!the directions are given in the cartesian basis
			 *
			 *       (0,1)
			 *         |
			 * (-1,0)--x--(1,0)
			 *         |
			 *      (0,-1)
			 *}*/
			this->dir_nn_[0](0) = 1.0;
			this->dir_nn_[0](1) = 0.0;

			this->dir_nn_[1](0) = 0.0;
			this->dir_nn_[1](1) = 1.0;

			this->dir_nn_[2](0) =-1.0;
			this->dir_nn_[2](1) = 0.0;

			this->dir_nn_[3](0) = 0.0;
			this->dir_nn_[3](1) =-1.0;

			this->x_[0] = this->dir_nn_[0]*0.5*(p_-q_)+this->dir_nn_[1]*0.5*(p_+q_);
			this->x_[0] = this->x_[0]/sqrt(this->x_[0].norm_squared())*0.001;

			Vector<double> x_loop(this->x_[0]);
			for(unsigned int i(1);i<this->n_;i++){
				this->x_[i] = this->x_[i-1] + this->dir_nn_[0];
				reset_pos_in_lattice(this->x_[i]);
				if(my::are_equal(this->x_[i],x_loop)){
					this->x_[i] += this->dir_nn_[1];
					reset_pos_in_lattice(this->x_[i]);
					x_loop = this->x_[i];
				}
				this->x_[i] = this->x_[i].chop();
			}

			this->equivalent_vertex_[0] = this->dir_nn_[0]*0.5*(p_-q_)+this->dir_nn_[1]*0.5*(p_+q_);
			this->equivalent_vertex_[1] = this->dir_nn_[1]*0.5*(p_-q_)+this->dir_nn_[2]*0.5*(p_+q_);
			this->equivalent_vertex_[2] = this->dir_nn_[3]*0.5*(p_-q_)+this->dir_nn_[0]*0.5*(p_+q_);

			if(this->unit_cell_allowed()){
				if(this->ref_(4)==2){ this->energy_obs(Vector<unsigned int>(1,2)); }
				else {
					this->ref_(4) = 0;
					this->status_ = 2;
				}

				/*!sets the bond energy if it has not been set yet*/
				if(this->obs_[0].nlinks() != this->J_.size()){
					if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
					else if (this->J_.size() == 2){
						Vector<double> tmp(this->J_);
						this->J_.set(this->obs_[0].nlinks());
						for (unsigned int i(0); i<this->J_.size();i++){
							if(this->obs_[0](i,3)){ this->J_(i) = tmp(1); } //vertical links
							else { this->J_(i) = tmp(0); } //horizontal links
						}

						/*!fix the names for the bond energy*/
						if(this->J_.size()==this->obs_[0].nlinks()){
							std::string tmp("Jp"+my::tostring(this->J_(1)));
							this->filename_.replace(this->filename_.find("Juniform"),8,tmp);
							this->path_.replace(this->path_.find("Juniform"),8,tmp);
						} else { std::cerr<<__PRETTY_FUNCTION__<<" : J_ has an incoherent size"<<std::endl; }
					} else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
				}
			}
		} else { std::cerr<<__PRETTY_FUNCTION__<<" required memory has not been allocated"<<std::endl; }
	}
}

template<typename Type>
void Square<Type>::draw_lattice(bool const& only_unit_cell, bool const& silent, bool const& only_lattice, Vector<double> const& uc_shift, std::string param, std::string title){
	Matrix<int> links(this->obs_[0].get_links());
	Vector<unsigned int> o(3,0);
	double max_bond_energy(0);
	for(unsigned int i(1);i<this->obs_.size();i++){
		switch(this->obs_[i].get_type()){
			case 1:
				{
					o(0)=i;
					for(unsigned int j(0);j<this->obs_[i].nval();j++){
						if(max_bond_energy < std::abs(this->obs_[i][j].get_x()/(this->m_*this->m_))){
							max_bond_energy = std::abs(this->obs_[i][j].get_x()/(this->m_*this->m_));
						}
					}
				}break;//bond energy
			case 2:{ o(1)=i; }break;//long range correlation
			case 3:{ o(2)=i; }break;//color occupation
		}
	}
	double max_chemical_potential(0);
	for(unsigned int j(0);j<this->n_;j++){
		if(max_chemical_potential < std::abs(my::real(this->H_(j,j)))){
			max_chemical_potential = std::abs(my::real(this->H_(j,j)));
		}
	}

	std::string arrow("-");
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	unsigned int s0;
	unsigned int s1;
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	Type t;
	double mu;
	double bond_energy;
	Vector<double> shift(2,0.0);
	Vector<unsigned int> loop(4);
	loop(0) = 0;
	loop(1) = 1;
	loop(2) = 2;
	loop(3) = 3;
	Matrix<double> uc(this->draw_unit_cell(uc_shift(0),uc_shift(1)));
	PSTricks ps(this->get_info_path(),this->filename_);
	ps.add("\\newcommand{\\wbg}[1]{\\setlength{\\fboxsep}{ 1pt}\\colorbox{white}{\\tiny{#1}}}");
	ps.begin(-20,-20,20,20,this->filename_);
	ps.polygon(uc,"linecolor=black,linestyle=dashed");
	if(only_unit_cell){
		shift(0) = uc(2,0)+1.0;
		for(unsigned int i(0);i<links.row();i++){
			s0 = links(i,0);
			s1 = links(i,1);
			xy0 = this->x_[s0];
			xy1 = this->x_[s1];

			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+this->dir_nn_[links(i,3)]).chop();
			} else { linestyle = "solid"; }

			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){
				/*Shows hopping amplitude, chemical potential and fluxes*/
				t = this->H_(s0,s1);
				if(std::abs(t)>1e-4){
					linewidth = my::tostring(std::abs(t))+"mm";
					if(my::real(t)>0){ color = "blue"; }
					else             { color = "red"; }
					if(my::are_equal(my::imag(t),0)){ arrow = "-"; }
					else {
						if(my::imag(-t)>0){ arrow = "->"; }
						else              { arrow = "<-"; }
					}

					if(links(i,4)){ linestyle = "dashed"; }
					else { linestyle = "solid"; }

					ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
					ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0, "\\wbg{"+my::tostring(my::round_nearest(std::abs(t),1000))+"}");
				}
				if(i%2){
					mu = my::real(this->H_(s0,s0));
					if(std::abs(mu)>1e-4){
						if(mu>0){ color = "cyan"; }
						else    { color = "magenta"; }
						ps.circle(xy0,sqrt(std::abs(mu)/max_chemical_potential*0.1),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
						ps.put(xy0(0),xy0(1), "\\wbg{"+my::tostring(my::round_nearest(std::abs(mu),1000))+"}");
					}
				} else { ps.put(xy0(0)+0.5,xy0(1)+0.5,this->flux_per_plaquette(s0,loop)); }

				/*Shows bond energy and color occupation*/
				xy0 += shift;
				xy1 += shift;
				if(o(0)){
					bond_energy = this->obs_[o(0)][links(i,2)].get_x()/(this->m_*this->m_);
					linewidth = my::tostring(std::abs(bond_energy))+"mm";
					if(std::abs(bond_energy)>1e-4){
						if(bond_energy>0){ color = "blue"; }
						else             { color = "red"; }
						ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					}
					ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0, "\\wbg{"+my::tostring(my::round_nearest(std::abs(bond_energy)/max_bond_energy,100))+"}");
				}
				if(i%2 && o(2)){
					Vector<double> p(this->N_);
					for(unsigned int j(0);j<this->N_;j++){ p(j) = this->obs_[o(2)][j+this->N_*links(i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
			}
		}
	} else {
		ps.polygon(this->cluster_vertex_,"linecolor=green");
		ps.linked_lines("-",this->draw_boundary(false),"linecolor=yellow");
		for(unsigned int i(0);i<links.row();i++){
			s0 = links(i,0);
			s1 = links(i,1);
			xy0 = this->x_[s0];
			xy1 = this->x_[s1];

			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+this->dir_nn_[links(i,3)]*1.2).chop();
				ps.put(xy1(0),xy1(1),"\\tiny{"+my::tostring(s1)+"}");
				xy1 = (xy0+this->dir_nn_[links(i,3)]).chop();
			} else { linestyle = "solid"; }

			/*Draws only the lattice, shows links and bc and indices*/
			linewidth = my::tostring(this->J_(i)/this->J_.max()*1)+"pt";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor=black,linestyle="+linestyle);
			if(i%2){
				ps.put(xy0(0)-0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
				ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\textcolor{green}{\\tiny{"+my::tostring(links(i,5))+"}}");
			}

			if(!only_lattice){
				/*Bond energy and color occupation*/
				shift = this->equivalent_vertex_[1]-this->equivalent_vertex_[2];
				xy0 += shift;
				xy1 += shift;
				if(o(0)){
					bond_energy = this->obs_[o(0)][links(i,2)].get_x()/(this->m_*this->m_);
					linewidth = my::tostring(std::abs(bond_energy))+"mm";
					if(std::abs(bond_energy)>1e-4){
						if(bond_energy>0){ color = "blue"; }
						else             { color = "red"; }
						ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					}
				}
				if(!(i%2) && o(2)){
					Vector<double> p(this->N_);
					for(unsigned int j(0);j<this->N_;j++){ p(j) = this->obs_[o(2)][j+this->N_*links(i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}

				/*Shows hopping amplitude, chemical potential and fluxes*/
				shift = this->equivalent_vertex_[0]-this->equivalent_vertex_[1];
				xy0 += shift;
				xy1 += shift;
				t = this->H_(s0,s1);
				if(std::abs(t)>1e-4){
					linewidth = my::tostring(std::abs(t))+"mm";
					if(my::real(t)>0){ color = "blue"; }
					else             { color = "red"; }
					if(my::are_equal(my::imag(t),0)){ arrow = "-"; }
					else {
						if(my::imag(-t)>0){ arrow = "->"; }
						else              { arrow = "<-"; }
					}
					ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
				}
				if(i%2){
					mu = my::real(this->H_(s0,s0))/M_PI;
					if(std::abs(mu)>1e-4){
						if(mu>0){ color = "cyan"; }
						else    { color = "magenta"; }
						ps.circle(xy0,sqrt(std::abs(mu)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
					}
				} else { ps.put(xy0(0)+0.5,xy0(1)+0.5,this->flux_per_plaquette(s0,loop)); }
			}
		}
	}
	/*Draws long range correlations over the lattice*/
	if(o(1)){ this->draw_long_range_correlations(ps,shift,this->obs_[o(1)]); }
	ps.end(silent,true,true);

	title = RST::math("J_p=")+my::tostring(this->J_(1)) + " : " + title;
	param += " -d:Jp " + my::tostring(this->J_(1));
	this->rst_file_set_default_info(param,title,"Jp"+my::tostring(this->J_(1))); 
}

template<typename Type>
std::string Square<Type>::extract_level_6(){
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->J_(1)<<" "<<this->obs_[0][0]<<" "<<this->ref_<<IOFiles::endl;

	return this->filename_;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Square<Type>::set_geometry(unsigned int const& n, unsigned int const& k, unsigned int const& spuc, unsigned int& ref3){
	if(n){
		p_ = sqrt(n);
		Matrix<double> tmp;
		if(!ref3 && my::are_equal(sqrt(n),p_)){
			tmp.set(5,2);
			q_ = 0;
		} else {
			for(unsigned int p(0);p<=sqrt(n);p++){
				for(unsigned int q(0);q<p+1;q++){
					if(p*p+q*q==n){
						tmp.set(5,2);
						p_ = p;
						q_ = q;
						p = n;
						q = n;
					}
				}
			}
			if(!q_){ ref3 = 0; }
		}
		if(tmp.ptr()){
			tmp(0,0) =-0.5*(p_-q_);
			tmp(0,1) =-0.5*(p_+q_);
			tmp(1,0) = 0.5*(p_+q_);
			tmp(1,1) =-0.5*(p_-q_);
			tmp(2,0) = 0.5*(p_-q_);
			tmp(2,1) = 0.5*(p_+q_);
			tmp(3,0) =-0.5*(p_+q_);
			tmp(3,1) = 0.5*(p_-q_);
			tmp(4,0) = tmp(0,0);
			tmp(4,1) = tmp(0,1);
			return tmp;
		}

		std::set<unsigned int> v;
		unsigned int m;
		for(unsigned int p(2);p<n;p++){
			for(unsigned int q(0);q<p+1;q++){
				m = p*p+q*q;
				if(!(m%spuc) && !(m%k) && m<2000){ v.insert(m); }
			}
		}
		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		for(auto const& n:v){ std::cerr<<"n="<<n<<std::endl; }
		std::cerr<<"n=p*p+q*q"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Square<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		double a(this->cluster_vertex_(0,0)*x(0)+this->cluster_vertex_(0,1)*x(1));
		double b(this->cluster_vertex_(1,0)*x(0)+this->cluster_vertex_(1,1)*x(1));
		if(a>0){
			if(b>0){ x += this->dir_nn_[1]*p_+this->dir_nn_[2]*q_; }
			else   { x += this->dir_nn_[0]*p_+this->dir_nn_[1]*q_; }
		} else {
			if(b>0){ x += this->dir_nn_[2]*p_+this->dir_nn_[3]*q_; }
			else   { x += this->dir_nn_[3]*p_+this->dir_nn_[0]*q_; }
		}
		reset_pos_in_lattice(x); //misplaced site (e.g. n=18, site 4)
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Square<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const {
	nn_dir = d;
	(void)(i);
	return this->dir_nn_[d];
}
/*}*/
#endif
