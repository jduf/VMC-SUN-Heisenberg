#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "System2DBis.hpp"

template<typename Type>
class Honeycomb: public System2DBis<Type>{
	public:
		/*!Constructor that organises the n=2L^2 sites (L integer)*/
		Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Honeycomb()=0;

	protected:
		void init_lattice();
		void set_obs(int nobs);

	private:
		double L_;

		Matrix<double> set_geometry(unsigned int const& n);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const;
};

/*{constructor*/
template<typename Type>
Honeycomb<Type>::Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2DBis<Type>(set_geometry((!this->obs_.size() || !this->obs_[0].nlinks())?this->n_:0),ab,spuc,3,filename)
{}

template<typename Type>
Honeycomb<Type>::~Honeycomb() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Honeycomb<Type>::init_lattice(){
	if(!this->obs_.size() || !this->obs_[0].nlinks()){
		/*{!the directions are given in the cartesian basis for the
		 * sublattice with even site number
		 *
		 * (-1,sqrt(3))/2
		 *       \
		 *        x--(1,0)
		 *       /
		 * (-1,-sqrt(3))/2
		 *
		 *  x = 0,2,4,...
		 *}*/
		this->dir_nn_[0](0) = 1.0;
		this->dir_nn_[0](1) = 0.0;

		this->dir_nn_[1](0) =-0.5;
		this->dir_nn_[1](1) = sqrt(3.0)/2.0;

		this->dir_nn_[2](0) =-0.5;
		this->dir_nn_[2](1) =-sqrt(3.0)/2.0;

		if(this->ref_(3)){ this->x_[0] = this->dir_nn_[2]; }
		else {
			this->x_[0] = this->dir_nn_[0]*(-0.5);
			this->x_[0](0)+= 0.01;
			this->x_[0](1)+= 0.01;
		}

		Vector<double> x_loop(this->x_[0]);
		bool check_if_loop(false);
		for(unsigned int i(1);i<this->n_;i++){
			if(this->ref_(3)){
				if(i%2){ this->x_[i] = this->x_[i-1] + this->dir_nn_[0]; }
				else   { this->x_[i] = this->x_[i-1] - this->dir_nn_[1]; }
			} else {
				if(i%2){ this->x_[i] = this->x_[i-1] + this->dir_nn_[1]; }
				else   { this->x_[i] = this->x_[i-1] - this->dir_nn_[2]; }
			}
			if(reset_pos_in_lattice(this->x_[i])){ check_if_loop = true; }
			if(check_if_loop && my::are_equal(this->x_[i],x_loop)){
				check_if_loop = false;
				if(this->ref_(3)){ this->x_[i](1) += sqrt(3.0); }
				else { this->x_[i] += this->dir_nn_[0]-this->dir_nn_[1]; }
				reset_pos_in_lattice(this->x_[i]);
				x_loop = this->x_[i];
			}
			this->x_[i] = this->x_[i].chop();
		}

		if(this->ref_(3)){
			this->boundary_vertex_[0] = (this->dir_nn_[3]+this->dir_nn_[4])*0.5 + (this->dir_nn_[3]+this->dir_nn_[4])*L_;
			this->boundary_vertex_[1] = (this->dir_nn_[3]+this->dir_nn_[4])*0.5 + (this->dir_nn_[0]+this->dir_nn_[5])*L_;
			this->boundary_vertex_[2] = (this->dir_nn_[3]+this->dir_nn_[4])*0.5 + (this->dir_nn_[0]+this->dir_nn_[1])*L_*2.0;
			this->boundary_vertex_[3] = (this->dir_nn_[3]+this->dir_nn_[4])*0.5 + (this->dir_nn_[1]+this->dir_nn_[2])*L_;
		} else {
			this->boundary_vertex_[0] = this->dir_nn_[3]*0.25 + this->dir_nn_[4]*L_;
			this->boundary_vertex_[1] = this->dir_nn_[3]*0.25 + this->dir_nn_[0]*L_;
			this->boundary_vertex_[2] = this->dir_nn_[3]*0.25 + this->dir_nn_[1]*L_*2.0;
			this->boundary_vertex_[3] = this->dir_nn_[3]*0.25 + this->dir_nn_[2]*L_;
		}

		Vector<unsigned int> l(2);
		l(0) = 3;
		l(1) = 0;
		this->set_nn_links(l);

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size()){
			if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
			else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
		}
	}
}

template<typename Type>
void Honeycomb<Type>::set_obs(int nobs){
	if(nobs<0){ nobs = 1; }
	if(nobs>0){/*bond energy (valid for Honeycomb0pp)*/
		unsigned int nlinks;
		unsigned int nval;
		nlinks = this->obs_[0].nlinks();
		nval = this->z_*this->spuc_/2;
		this->obs_.push_back(Observable("Bond energy",1,nval,nlinks));
		this->obs_[1].remove_links();
		for(unsigned int i(0);i<nlinks;i++){
			this->obs_[0](i,2) = this->obs_[0](i,5)/2*3;
			this->obs_[0](i,2)+=(this->obs_[0](i,6)-1)/2;
		}
	}
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Honeycomb<Type>::set_geometry(unsigned int const& n){
	if(n){
		L_ = sqrt(n/2.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2);
			Matrix<double> tmp(7,2);
			tmp(0,0) =-0.5*L_;
			tmp(0,1) =-a*L_;
			tmp(1,0) =-tmp(0,0);
			tmp(1,1) = tmp(0,1);
			tmp(2,0) = L_;
			tmp(2,1) = 0;
			tmp(3,0) =-tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(2,0);
			tmp(5,1) = tmp(2,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}
		L_ = sqrt(n/6.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2.0);
			Matrix<double> tmp(7,2);
			tmp(0,0) = 0;
			tmp(0,1) =-2.0*a*L_;
			tmp(1,0) = 1.5*L_;
			tmp(1,1) =-a*L_;
			tmp(2,0) = tmp(1,0);
			tmp(2,1) =-tmp(1,1);
			tmp(3,0) = tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(1,0);
			tmp(5,1) = tmp(1,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}
		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		for(unsigned int l(2);l<10;l++){ std::cerr<<"n="<<2*l*l<<" or "<<6*l*l<<std::endl; }
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
				if(std::abs(t)>1){ x-= this->dir_nn_[0]*L_*3.0; }
				else {
					if(t>0){       x+= this->dir_nn_[2]*L_*3.0; }
					else   {       x+= this->dir_nn_[1]*L_*3.0; }
				}
			} else {
				if(std::abs(t)>1){ x+= this->dir_nn_[0]*L_*3.0; }
				else {
					if(t>0){       x-= this->dir_nn_[2]*L_*3.0; }
					else   {       x-= this->dir_nn_[1]*L_*3.0; }
				}
			}
		} else {
			double t(tan(M_PI/3.0)*x(0)/x(1));
			if(x(1)>0){
				if(std::abs(t)<1){ x-=(this->dir_nn_[1]-this->dir_nn_[2])*L_; }
				else {
					if(t>0){       x-=(this->dir_nn_[0]-this->dir_nn_[2])*L_; }
					else   {       x-=(this->dir_nn_[1]-this->dir_nn_[0])*L_; }
				}
			} else {
				if(std::abs(t)<1){ x+=(this->dir_nn_[1]-this->dir_nn_[2])*L_; }
				else {
					if(t>0){       x+=(this->dir_nn_[0]-this->dir_nn_[2])*L_; }
					else   {       x+=(this->dir_nn_[1]-this->dir_nn_[0])*L_; }
				}
			}
		}
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Honeycomb<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const {
	if(i%2){ return -this->dir_nn_[d]; }
	else { return this->dir_nn_[d]; }
}
/*}*/
#endif
