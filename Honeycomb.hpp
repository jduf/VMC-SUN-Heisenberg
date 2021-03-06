#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "System2D.hpp"

template<typename Type>
class Honeycomb: public System2D<Type>{
	public:
		/*!Constructor that organises the n=2L^2 or 6L^2 sites (L integer)*/
		Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Honeycomb()=0;

	protected:
		void init_lattice();
		/*Draw the lattice inside a PSTricks file*/
		void draw_lattice(bool const& only_unit_cell, bool const& silent, bool const& only_lattice, Vector<double> const& uc_shift, std::string const& param, std::string const& title);
		/*!Create the long range correlation observables*/
		void long_range_correlations_obs();

		virtual void param_fit_therm_limit(std::string& f, std::string& param, std::string& range);
		virtual std::string extract_level_2();

	private:
		double L_;

		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& k, unsigned int const& spuc);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const;
};

/*{constructor*/
template<typename Type>
Honeycomb<Type>::Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry((this->ref_(4)?this->n_:0),this->N_/this->m_,spuc),ab,spuc,3,6,filename)
{}

template<typename Type>
Honeycomb<Type>::~Honeycomb() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Honeycomb<Type>::init_lattice(){
	if(!this->ref_(4)){ this->status_ = 2; }
	else {
		if(this->dir_nn_){
			/*{!The directions are given in the cartesian basis
			 *
			 * (-1,sqrt(3))/2  (1,sqrt(3))/2
			 *               \ /
			 *        (-1,0)--x--(1,0)
			 *               / \
			 * (-1,-sqrt(3))/2 (1,-sqrt(3))/2
			 *}*/
			this->dir_nn_[0](0) = 1.0;
			this->dir_nn_[0](1) = 0.0;

			this->dir_nn_[1](0) = 0.5;
			this->dir_nn_[1](1) = sqrt(3.0)/2.0;

			this->dir_nn_[2](0) =-0.5;
			this->dir_nn_[2](1) = sqrt(3.0)/2.0;

			this->dir_nn_[3](0) =-1.0;
			this->dir_nn_[3](1) = 0.0;

			this->dir_nn_[4](0) =-0.5;
			this->dir_nn_[4](1) =-sqrt(3.0)/2.0;

			this->dir_nn_[5](0) = 0.5;
			this->dir_nn_[5](1) =-sqrt(3.0)/2.0;

			if(this->ref_(3)){ this->x_[0] = this->dir_nn_[4]; }
			else {
				this->x_[0] = this->dir_nn_[3]*0.5;
				this->x_[0](0)+= 0.01;
				this->x_[0](1)+= 0.01;
			}

			Vector<double> x_loop(this->x_[0]);
			for(unsigned int i(1);i<this->n_;i++){
				if(this->ref_(3)){
					if(i%2){ this->x_[i] = this->x_[i-1] + this->dir_nn_[0]; }
					else   { this->x_[i] = this->x_[i-1] + this->dir_nn_[5]; }
				} else {
					if(i%2){ this->x_[i] = this->x_[i-1] + this->dir_nn_[2]; }
					else   { this->x_[i] = this->x_[i-1] + this->dir_nn_[1]; }
				}
				reset_pos_in_lattice(this->x_[i]);
				if(my::are_equal(this->x_[i],x_loop)){
					if(this->ref_(3)){ this->x_[i](1) += sqrt(3.0); }
					else { this->x_[i] += this->dir_nn_[0]+this->dir_nn_[5]; }
					reset_pos_in_lattice(this->x_[i]);
					x_loop = this->x_[i];
				}
				this->x_[i] = this->x_[i].chop();
			}

			if(this->ref_(3)){
				this->equivalent_vertex_[0] = (this->dir_nn_[4]+this->dir_nn_[3])*L_ + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
				this->equivalent_vertex_[1] = (this->dir_nn_[0]+this->dir_nn_[5])*L_ + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
				this->equivalent_vertex_[2] = (this->dir_nn_[2]+this->dir_nn_[1])*L_ + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
			} else {
				this->equivalent_vertex_[0] = this->dir_nn_[4]*L_ + (this->dir_nn_[2] + this->dir_nn_[1])/2.0;
				this->equivalent_vertex_[1] = this->dir_nn_[0]*L_ + (this->dir_nn_[2] + this->dir_nn_[1])/2.0;
				this->equivalent_vertex_[2] = this->dir_nn_[2]*L_ + (this->dir_nn_[2] + this->dir_nn_[1])/2.0;
			}

			if(this->unit_cell_allowed()){
				if(this->ref_(4)==2){
					Vector<unsigned int> l(2);
					l(0) = 3;
					l(1) = 0;
					this->energy_obs(l);
				} else {
					this->ref_(4) = 0;
					this->status_ = 2;
				}

				/*!sets the bond energy if it has not been set yet*/
				if(this->obs_[0].nlinks() != this->J_.size()){
					if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
				}
			}
		} else { std::cerr<<__PRETTY_FUNCTION__<<" required memory has not been allocated"<<std::endl; }
	}
}

template<typename Type>
void Honeycomb<Type>::long_range_correlations_obs(){
	Vector<double>* dx(new Vector<double>[this->n_]);
	for(unsigned int i(0);i<this->n_;i++){ dx[i] = this->x_[i]-this->x_[0]; }

	this->obs_.push_back(Observable("Long range correlations",2,this->n_,this->n_*this->n_/2));
	for(unsigned int i(0);i<this->n_/2;i++){
		for(unsigned int j(0);j<this->n_;j++){
			this->obs_.back()(i*this->n_+j,0) = 2*i;
			this->obs_.back()(i*this->n_+j,1) = this->site_index(this->x_[2*i]+dx[j]);
			this->obs_.back()(i*this->n_+j,2) = j;
		}
	}
	delete[] dx;
}

template<typename Type>
void Honeycomb<Type>::draw_lattice(bool const& only_unit_cell, bool const& silent, bool const& only_lattice, Vector<double> const& uc_shift, std::string const& param, std::string const& title){
	Matrix<int> links(this->obs_[0].get_links());
	Vector<unsigned int> o(3,0);
	for(unsigned int i(1);i<this->obs_.size();i++){
		switch(this->obs_[i].get_type()){
			case 1:{ o(0)=i; }break;//bond energy
			case 2:{ o(1)=i; }break;//long range correlation
			case 3:{ o(2)=i; }break;//color occupation
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
	Vector<unsigned int> loop(6);
	loop(0) = 0;
	loop(1) = 1;
	loop(2) = 2;
	loop(3) = 3;
	loop(4) = 4;
	loop(5) = 5;
	Matrix<double> uc(this->draw_unit_cell(uc_shift(0),uc_shift(1)));
	PSTricks ps(this->get_info_path(),this->filename_);
	ps.add("\\newcommand{\\wbg}[1]{\\setlength{\\fboxsep}{ 1pt}\\colorbox{white}{\\tiny{#1}}}");
	ps.begin(-20,-20,20,20,this->filename_);
	ps.polygon(uc,"linecolor=black,linestyle=dashed");
	if(only_unit_cell){
		shift(0) = uc(1,0)-uc(0,0)+1.5;
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
					ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0, "\\wbg{"+my::tostring(my::round_nearest(std::abs(t),1000))+"}");
				}
				mu = my::real(this->H_(s0,s0));
				if(std::abs(mu)>1e-4){
					std::cout<<__PRETTY_FUNCTION__<<" : the display of the chemical potential might be wrong, need to check the code"<<std::endl;
					if(mu>0){ color = "cyan"; }
					else    { color = "magenta"; }
					ps.circle(xy0,sqrt(std::abs(mu)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
				}
				if(!(i%3)){ ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.9,this->flux_per_plaquette(s0,loop)); }

				/*Shows bond energy and color occupation*/
				xy0 += shift;
				xy1 += shift;
				if(o(0)){
					bond_energy = this->obs_[o(0)][links(i,2)].get_x()/(this->m_*this->m_*this->J_(i));
					linewidth = my::tostring(std::abs(bond_energy))+"mm";
					if(std::abs(bond_energy)>1e-4){
						if(bond_energy>0){ color = "blue"; }
						else             { color = "red"; }
						ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					}
					ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0, "\\wbg{"+my::tostring(my::round_nearest(std::abs(bond_energy),100))+"}");
				}
				if(i%2 && o(2)){
					std::cout<<__PRETTY_FUNCTION__<<" : the display of the color occupation might be wrong, need to check the code"<<std::endl;
					Vector<double> p(this->N_);
					for(unsigned int j(0);j<this->N_;j++){ p(j) = this->obs_[o(2)][j+this->N_*links(i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
			}
		}
	} else {
		ps.polygon(this->cluster_vertex_,"linecolor=green");
		ps.linked_lines("-",this->draw_boundary(false),"linecolor=yellow");
		Vector<double> xy0tmp;
		Vector<double> xy1tmp;
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
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth=1pt,linecolor=black,linestyle="+linestyle);
			if(!(i%3)){
				xy0tmp = xy0;
				xy1tmp = xy1;
				xy0tmp -=  this->dir_nn_[this->obs_[0](i,3)]*0.2;
				xy1tmp -=  this->dir_nn_[this->obs_[0](i,3)]*0.2;
				ps.put(xy0tmp(0),xy0tmp(1),"\\tiny{"+my::tostring(this->obs_[0](i,0))+"}");
				xy0tmp -=  this->dir_nn_[this->obs_[0](i,3)]*0.2;
				xy1tmp -=  this->dir_nn_[this->obs_[0](i,3)]*0.2;
				ps.put(xy0tmp(0),xy0tmp(1),"\\textcolor{green}{\\tiny{"+my::tostring(this->obs_[0](i,5))+"}}");

				xy0tmp +=  this->dir_nn_[this->obs_[0](i,3)]*0.6;
				xy1tmp +=  this->dir_nn_[this->obs_[0](i,3)]*0.6;
				ps.put(xy1tmp(0),xy1tmp(1),"\\tiny{"+my::tostring(this->obs_[0](i,1))+"}");
				xy0tmp +=  this->dir_nn_[this->obs_[0](i,3)]*0.2;
				xy1tmp +=  this->dir_nn_[this->obs_[0](i,3)]*0.2;
				ps.put(xy1tmp(0),xy1tmp(1),"\\textcolor{green}{\\tiny{"+my::tostring(this->obs_[0](i,6))+"}}");
			}

			if(!only_lattice){
				/*Bond energy and color occupation*/
				shift = this->equivalent_vertex_[1]-this->equivalent_vertex_[2];
				xy0 += shift;
				xy1 += shift;
				if(o(0)){
					bond_energy = this->obs_[o(0)][links(i,2)].get_x()/(this->m_*this->m_*this->J_(i));
					linewidth = my::tostring(std::abs(bond_energy))+"mm";
					if(std::abs(bond_energy)>1e-4){
						if(bond_energy>0){ color = "blue"; }
						else             { color = "red"; }
						ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					}
				}
				if(i%2 && o(2)){
					std::cout<<__PRETTY_FUNCTION__<<" : the display of the color occupation might be wrong, need to check the code"<<std::endl;
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
				mu = my::real(this->H_(s0,s0));
				if(std::abs(mu)>1e-4){
					std::cout<<__PRETTY_FUNCTION__<<" : the display of the color occupation might be wrong, need to check the code"<<std::endl;
					if(mu>0){ color = "cyan"; }
					else    { color = "magenta"; }
					ps.circle(xy0,sqrt(std::abs(mu)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
				}
				if(!(i%3)){ ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.9,this->flux_per_plaquette(s0,loop)); }
			}
		}
	}
	/*draws long range correlations over the lattice*/
	if(o(1)){ this->draw_long_range_correlations(ps,shift,this->obs_[o(1)]); }
	ps.end(silent,true,true);

	this->rst_file_set_default_info(param,title);
}

template<typename Type>
void Honeycomb<Type>::param_fit_therm_limit(std::string& f, std::string& param, std::string& range){
	f="f(x) = a+b*x*x";
	param = "a,b";
	range = "[0:0.025]";
}

template<typename Type>
std::string Honeycomb<Type>::extract_level_2(){
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_);
	gp+="f(x) = a+b*x*x";
	gp+="g(x,a,b) = a+b*x*x";
	gp+="set fit quiet";
	gp+="fit [0:0.025] f(x) '"+this->filename_+".dat' u ($1<0 && $12==4 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) yerror via a,b";
	gp+="a0=a";
	gp+="b0=b";
	gp+="fit [0:0.025] f(x) '"+this->filename_+".dat' u ($1>0 && $12==4 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) yerror via a,b";
	gp+="a1=a";
	gp+="b1=b";
	gp+="fit [0:0.025] f(x) '"+this->filename_+".dat' u ($1<0 && $12==3 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) yerror via a,b";
	gp+="a2=a";
	gp+="b2=b";
	gp+="fit [0:0.025] f(x) '"+this->filename_+".dat' u ($1>0 && $12==3 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) yerror via a,b";
	gp+="a3=a";
	gp+="b3=b";
	gp+="set print \'../"+this->sim_.substr(0,this->sim_.size()-1)+".dat\' append";
	gp+="print " + my::tostring(this->N_) + "," + my::tostring(this->m_) + ",a0,a1,a2,a3";
	gp.range("x","0","");
	gp.label("x","$\\frac{ 1}{n}$");
	gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
	gp+="plot '"+this->filename_+".dat' u ($1<0 && $12==4 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) w e lc 7 notitle,\\";
	gp+="     '"+this->filename_+".dat' u ($1>0 && $12==4 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) w e lc 4 notitle,\\";
	gp+="     '"+this->filename_+".dat' u ($1<0 && $12==3 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) w e lc 3 notitle,\\";
	gp+="     '"+this->filename_+".dat' u ($1>0 && $12==3 ?1.0/$4:1/0):($6/($2*$2)):($7/($2*$2)) w e lc 2 notitle,\\";
	gp+="     [0:0.025] g(x,a0,b0) lc 7 t sprintf('$\\pi\\pi\\pi,\\,%f$',a0),\\";
	gp+="     [0:0.025] g(x,a1,b1) lc 4 t sprintf('$\\pi00,\\,%f$',a1),\\";
	gp+="     [0:0.025] g(x,a2,b2) lc 3 t sprintf('$000,\\,%f$',a2),\\";
	gp+="     [0:0.025] g(x,a3,b3) lc 2 t sprintf('$0\\pi\\pi,\\,%f$',a3)";
	gp.save_file();
	gp.create_image(true,"png");

	return this->filename_;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Honeycomb<Type>::set_geometry(unsigned int const& n, unsigned int const& k, unsigned int const& spuc){
	if(n){
		L_ = sqrt(n/2.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2.0);
			Matrix<double> tmp(7,2);
			tmp(0,0) =-0.5*L_;
			tmp(0,1) =-a*L_;
			tmp(1,0) =-tmp(0,0);
			tmp(1,1) = tmp(0,1);
			tmp(2,0) = L_;
			tmp(2,1) = 0.0;
			tmp(3,0) =-tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(2,0);
			tmp(5,1) =-tmp(2,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}
		L_ = sqrt(n/6.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2.0);
			Matrix<double> tmp(7,2);
			tmp(0,0) = 0.0;
			tmp(0,1) =-2.0*a*L_;
			tmp(1,0) = 1.5*L_;
			tmp(1,1) =-a*L_;
			tmp(2,0) = tmp(1,0);
			tmp(2,1) =-tmp(1,1);
			tmp(3,0) =-tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(1,0);
			tmp(5,1) = tmp(1,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}

		std::set<unsigned int> v;
		unsigned int m;
		for(unsigned int i(2);i<n;i++){
			m = 6*i*i;
			if(!(m%spuc) && !(m%k) && m<2000){ v.insert(m); }
			m = 2*i*i;
			if(!(m%spuc) && !(m%k) && m<2000){ v.insert(m); }
		}
		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		for(auto const& n:v){ std::cerr<<"n="<<n<<std::endl; }
		std::cerr<<"n=2*l*l or 6*l*l"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Honeycomb<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(this->ref_(3)){
			double t(tan(M_PI/6.0)*x(0)/x(1));
			if(x(0)>0){
				if(std::abs(t)>1){ x+= this->dir_nn_[3]*L_*3.0; }
				else {
					if(t>0){       x+= this->dir_nn_[4]*L_*3.0; }
					else   {       x+= this->dir_nn_[2]*L_*3.0; }
				}
			} else {
				if(std::abs(t)>1){ x+= this->dir_nn_[0]*L_*3.0; }
				else {
					if(t>0){       x+= this->dir_nn_[1]*L_*3.0; }
					else   {       x+= this->dir_nn_[5]*L_*3.0; }
				}
			}
		} else {
			double t(tan(M_PI/3.0)*x(0)/x(1));
			if(x(1)>0){
				if(std::abs(t)<1){ x+= (this->dir_nn_[4]+this->dir_nn_[5])*L_; }
				else {
					if(t>0){       x+= (this->dir_nn_[3]+this->dir_nn_[4])*L_; }
					else   {       x+= (this->dir_nn_[5]+this->dir_nn_[0])*L_; }
				}
			} else {
				if(std::abs(t)<1){ x+= (this->dir_nn_[2]+this->dir_nn_[1])*L_; }
				else {
					if(t>0){       x+= (this->dir_nn_[0]+this->dir_nn_[1])*L_; }
					else   {       x+= (this->dir_nn_[2]+this->dir_nn_[3])*L_; }
				}
			}
		}
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Honeycomb<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const {
	nn_dir = 2*d+(i%2);
	return this->dir_nn_[nn_dir];
}
/*}*/
#endif
