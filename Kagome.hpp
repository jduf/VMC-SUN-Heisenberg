#ifndef DEF_KAGOME
#define DEF_KAGOME

#include "System2D.hpp"

template<typename Type>
class Kagome: public System2D<Type>{
	public:
		/*!Constructor that organises the n=9L^2 or 9L^2 sites (L integer)*/
		Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Kagome()=0;

	protected:
		void init_lattice();
		void draw_lattice();
		void compute_long_range_correlation();

	private:
		double L_;

		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& spuc);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const;
};

/*{constructor*/
template<typename Type>
Kagome<Type>::Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry((this->ref_(4)?this->n_:0),spuc),ab,spuc,4,6,filename)
{}

template<typename Type>
Kagome<Type>::~Kagome() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Kagome<Type>::init_lattice(){
	if(!this->ref_(4)){ this->status_ = 2; }
	else {
		if(this->dir_nn_){
			/*{!the directions are given in the cartesian basis
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

			this->x_[0].set(2);
			this->x_[0](0) =-0.5001;
			if(this->ref_(3)){ this->x_[0](1) =-0.5001; }
			else             { this->x_[0](1) =-0.0001; }

			Vector<double> x_loop(this->x_[0]);
			for(unsigned int i(1);i<this->n_;i++){
				switch(i%3){
					case 0: { this->x_[i] = this->x_[i-1] + this->dir_nn_[0]; }break;
					case 1: { this->x_[i] = this->x_[i-1] + this->dir_nn_[1]; }break;
					case 2: { this->x_[i] = this->x_[i-1] + this->dir_nn_[5]; }break;
				}
				reset_pos_in_lattice(this->x_[i]);
				if(my::are_equal(this->x_[i],x_loop)){
					this->x_[i] += this->dir_nn_[1]*2.0;
					reset_pos_in_lattice(this->x_[i]);
					x_loop = this->x_[i];
				}
				this->x_[i] = this->x_[i].chop();
			}

			if(this->ref_(3)){
				this->equivalent_vertex_[0] = this->dir_nn_[4]*L_*2.0 + this->dir_nn_[3]*0.2;
				this->equivalent_vertex_[1] = this->dir_nn_[0]*L_*2.0 + this->dir_nn_[3]*0.2;
				this->equivalent_vertex_[2] = this->dir_nn_[2]*L_*2.0 + this->dir_nn_[3]*0.2;
			} else {
				this->equivalent_vertex_[0] = (this->dir_nn_[3]+this->dir_nn_[4])*L_/1.5 + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
				this->equivalent_vertex_[1] = (this->dir_nn_[5]+this->dir_nn_[0])*L_/1.5 + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
				this->equivalent_vertex_[2] = (this->dir_nn_[1]+this->dir_nn_[2])*L_/1.5 + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
			}

			if(this->unit_cell_allowed()){
				this->status_ = 2;

				if(this->ref_(4)==2){ this->create_energy_obs(Vector<unsigned int>(3,2)); }
				else { this->ref_(4) = 0; }

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
void Kagome<Type>::draw_lattice(){
	Matrix<int> links(this->obs_[0].get_links());
	Vector<unsigned int> o(3,0);
	for(unsigned int i(1);i<this->obs_.size();i++){
		switch(this->obs_[i].get_type()){
			case 1:{ o(0)=i; }break;//bond energy
			case 2:{ o(1)=i; }break;//long range correlation
			case 3:{ o(2)=i; }break;//color occupation
		}
	}

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(this->get_info_path(),this->filename_);
	ps.begin(-20,-20,20,20,this->filename_);
	ps.polygon(this->cluster_vertex_,"linecolor=green");
	Matrix<double> uc(this->draw_unit_cell(-1.0,-sqrt(3.0)/4.0));
	ps.polygon(uc,"linecolor=black");
	ps.linked_lines("-",this->draw_boundary(false),"linecolor=yellow");

	unsigned int s0;
	unsigned int s1;
	/*draws only the lattice, shows links and bc*/
	for(unsigned int i(0);i<links.row();i++){
		s0 = links(i,0);
		xy0 = this->x_[s0];
		s1 = links(i,1);
		xy1 = this->x_[s1];

		if((xy0-xy1).norm_squared()>1.0001){
			linestyle = "dashed";
			xy1 = (xy0+this->dir_nn_[links(i,3)]*1.2).chop();
			ps.put(xy1(0),xy1(1),"\\tiny{"+my::tostring(s1)+"}");
			xy1 = (xy0+this->dir_nn_[links(i,3)]).chop();
		} else { linestyle = "solid"; }

		ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth=1pt,linecolor=black,linestyle="+linestyle);
		if(i%2){
			switch(links(i,5)%3){
				case 0:
					{
						ps.put(xy0(0)-0.2,xy0(1)+0.2,"\\tiny{"+my::tostring(s0)+"}");
						ps.put(xy0(0)+0.2,xy0(1)-0.2,"\\textcolor{green}{\\tiny{"+my::tostring(links(i,5))+"}}"); 
					}break;
				case 1: 
					{
						ps.put(xy0(0)-0.2,xy0(1),"\\tiny{"+my::tostring(s0)+"}"); 
						ps.put(xy0(0)+0.2,xy0(1),"\\textcolor{green}{\\tiny{"+my::tostring(links(i,5))+"}}"); 
					}break;
				case 2: 
					{
						ps.put(xy0(0)-0.2,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); 
						ps.put(xy0(0)+0.2,xy0(1)+0.2,"\\textcolor{green}{\\tiny{"+my::tostring(links(i,5))+"}}"); 
					}break;
			}
		}
	}
	/*draws long range correlations over the lattice*/
	if(o(1)){ this->draw_long_range_correlation(ps,this->obs_[o(1)]); }

	Vector<double> shift(2,0.0);
	if(o(0) || o(2)){
		/*unit cell, shows bond energy and color occupation*/
		double be;
		//Vector<double> shift(equivalent_vertex_[0]+equivalent_vertex_[1]);
		//ps.polygon(draw_unit_cell(shift(0)+0.5,shift(1)+0.5),"linecolor=black");
		for(unsigned int i(0);i<links.row();i++){
			xy0 = this->x_[links(i,0)];
			xy1 = this->x_[links(i,1)];

			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1)))
			{
				xy0 += shift;
				xy1 += shift;
				if(o(0)){
					be = this->obs_[o(0)][links(i,2)].get_x();
					linewidth = my::tostring(std::abs(be))+"mm";
					if(std::abs(be)>1e-4){
						if(be>0){ color = "blue"; }
						else    { color = "red"; }
						ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					}
				}
				if(i%2 && o(2)){
					Vector<double> p(this->N_);
					for(unsigned int j(0);j<this->N_;j++){ p(j) = this->obs_[o(2)][j+this->N_*links(i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
			}
		}
	}

	/*unit cell, shows hopping amplitude, chemical potential and fluxes*/
	Type t;
	double mu;
	double flux;
	double sign;
	unsigned long long a;
	unsigned long long b;
	std::string arrow("-");
	//shift = equivalent_vertex_[0]+equivalent_vertex_[2];
	//ps.polygon(draw_unit_cell(shift(0)+0.5,shift(1)+0.5),"linecolor=black");
	for(unsigned int i(0);i<links.row();i++){
		s0 = links(i,0);
		xy0 = this->x_[s0];
		s1 = links(i,1);
		xy1 = this->x_[s1];

		if((xy0-xy1).norm_squared()>1.0001){
			linestyle = "dashed";
			xy1 = (xy0+this->dir_nn_[links(i,3)]*1.2).chop();
			ps.put(xy1(0),xy1(1),"\\tiny{"+my::tostring(s1)+"}");
			xy1 = (xy0+this->dir_nn_[links(i,3)]).chop();
		} else { linestyle = "solid"; }

		if(!my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1)))
		{
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
				ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			mu = my::real(this->H_(s0,s0));
			if(std::abs(mu)>1e-4){
				if(mu>0){ color = "cyan"; }
				else    { color = "magenta"; }
				ps.circle(xy0,sqrt(std::abs(mu)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}

			switch(links(i,5)%3){
				case 0:
					{
						if(links(i,3)==0){
							flux = 0.0;
							unsigned int j(0);
							do {
								xy0 += this->dir_nn_[2*j];
								s1 = this->site_index(xy0);
								flux += std::arg(-this->H_(s0,s1));
								s0 = s1;
							} while (++j<3);
							flux /= M_PI;
							if(!my::are_equal(flux,0.0,this->eq_prec_,this->eq_prec_)){
								if(my::to_fraction(flux,a,b,sign) && b!=1){
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
								} else {
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(my::chop(flux))+"}");
								}
							}
						}
					}break;
				case 1:
					{
						if(links(i,3)==1){
							flux = 0.0;
							unsigned int j(0);
							do {
								xy0 += this->dir_nn_[2*j+1];
								s1 = this->site_index(xy0);
								flux += std::arg(-this->H_(s0,s1));
								s0 = s1;
							} while (++j<3);
							flux /= M_PI;
							if(!my::are_equal(flux,0.0,this->eq_prec_,this->eq_prec_)){
								if(my::to_fraction(flux,a,b,sign) && b!=1){
									ps.put(xy0(0),(xy0(1)+xy1(1))/2.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
								} else {
									ps.put(xy0(0),(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(my::chop(flux))+"}");
								}
							}
						}
					}break;
				case 2:
					{
						if(links(i,3)==0){
							flux = 0.0;
							unsigned int j(0);
							do {
								xy0 += this->dir_nn_[j];
								s1 = this->site_index(xy0);
								flux += std::arg(-this->H_(s0,s1));
								s0 = s1;
							} while (++j<6);
							flux /= M_PI;
							if(!my::are_equal(flux,0.0,this->eq_prec_,this->eq_prec_)){
								if(my::to_fraction(flux,a,b,sign) && b!=1){
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+1.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
								} else {
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+1.0,"\\tiny{"+my::tostring(my::chop(flux))+"}");
								}
							}
						}
					}break;
			}
		}
	}
	ps.end(true,true,true);
}

template<typename Type>
void Kagome<Type>::compute_long_range_correlation(){
	Vector<double>* dx(new Vector<double>[this->n_]);
	for(unsigned int i(0);i<this->n_;i++){ dx[i] = this->x_[i]-this->x_[0]; }

	this->obs_.push_back(Observable("Long range correlations",2,this->n_,this->n_*this->n_/3));
	for(unsigned int i(0);i<this->n_/3;i++){
		for(unsigned int j(0);j<this->n_;j++){
			this->obs_.back()(i*this->n_+j,0) = 3*i;
			this->obs_.back()(i*this->n_+j,1) = this->site_index(this->x_[3*i]+dx[j]);
			this->obs_.back()(i*this->n_+j,2) = j;
		}
	}
	delete[] dx;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Kagome<Type>::set_geometry(unsigned int const& n, unsigned int const& spuc){
	if(n){
		L_ = sqrt(n/9.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2);
			Matrix<double> tmp(7,2);
			tmp(0,0) =-L_;
			tmp(0,1) =-L_*a*2.0;
			tmp(1,0) =-tmp(0,0);
			tmp(1,1) = tmp(0,1);
			tmp(2,0) = L_*2.0;
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
		L_ = sqrt(n/3.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2.0);
			Matrix<double> tmp(7,2);
			tmp(0,0) = 0.0;
			tmp(0,1) =-L_/a;
			tmp(1,0) = L_;
			tmp(1,1) = tmp(0,1)/2.0;
			tmp(2,0) = tmp(1,0);
			tmp(2,1) =-tmp(1,1);
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

		std::set<unsigned int> v;
		unsigned int m;
		for(unsigned int i(2);i<20;i++){
			m = 3*i*i;
			if(!(m%spuc)){ v.insert(m); }
			m = (3*i)*(3*i);
			if(!(m%spuc)){ v.insert(m); }
		}
		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		for(auto const& n:v){ std::cerr<<"n="<<n<<std::endl; }
		std::cerr<<"n=9*l*l or 3*l*l"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Kagome<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(this->ref_(3)){
			double t(tan(M_PI/3.0)*x(0)/x(1));
			if(x(1)>0){
				if(std::abs(t)<1){ x+= (this->dir_nn_[4]+this->dir_nn_[5])*L_*2.0; }
				else {
					if(t>0){       x+= (this->dir_nn_[4]+this->dir_nn_[3])*L_*2.0; }
					else   {       x+= (this->dir_nn_[5]+this->dir_nn_[0])*L_*2.0; }
				}
			} else {
				if(std::abs(t)<1){ x+= (this->dir_nn_[1]+this->dir_nn_[2])*L_*2.0; }
				else {
					if(t>0){       x+= (this->dir_nn_[1]+this->dir_nn_[0])*L_*2.0; }
					else   {       x+= (this->dir_nn_[2]+this->dir_nn_[3])*L_*2.0; }
				}
			}
		} else {
			double t(tan(M_PI/6.0)*x(0)/x(1));
			if(x(0)>0){
				if(std::abs(t)>1){ x+= this->dir_nn_[3]*L_*2.0; }
				else {
					if(t>0){       x+= this->dir_nn_[4]*L_*2.0; }
					else   {       x+= this->dir_nn_[2]*L_*2.0; }
				}
			} else {
				if(std::abs(t)>1){ x+= this->dir_nn_[0]*L_*2.0; }
				else {
					if(t>0){       x+= this->dir_nn_[1]*L_*2.0; }
					else   {       x+= this->dir_nn_[5]*L_*2.0; }
				}
			}
		}
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Kagome<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const {
	switch(i%3){
		case 0:
			{
				switch(d){
					case 0: { nn_dir = 0; }break;
					case 1: { nn_dir = 1; }break;
					case 2: { nn_dir = 3; }break;
					case 3: { nn_dir = 4; }break;
				}
			}break;
		case 1:
			{
				switch(d){
					case 0: { nn_dir = 1; }break;
					case 1: { nn_dir = 2; }break;
					case 2: { nn_dir = 4; }break;
					case 3: { nn_dir = 5; }break;
				}
			}break;
		case 2:
			{
				switch(d){
					case 0: { nn_dir = 0; }break;
					case 1: { nn_dir = 2; }break;
					case 2: { nn_dir = 3; }break;
					case 3: { nn_dir = 5; }break;
				}
			}break;
	}
	return this->dir_nn_[nn_dir];
}
/*}*/
#endif
