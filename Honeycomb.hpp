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
		void set_obs(int nobs);

	private:
		double L_;
		unsigned int lattice_type_;

		Matrix<double> set_geometry(unsigned int const& n);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const;
};

/*{constructor*/
template<typename Type>
Honeycomb<Type>::Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2DBis<Type>(set_geometry(this->n_),ab,spuc,3,filename)
{
	if(this->status_==2){
		/*{!the directions are given for the sublattice with even site number
		 * in the cartesian basis
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

		if(lattice_type_){ this->x_[0] = this->dir_nn_[2]; }
		else {
			this->x_[0] = this->dir_nn_[0]*(-0.5);
			this->x_[0](0)+= 0.01;
			this->x_[0](1)+= 0.01;
		}

		//PSTricks ps("./","test");
		//ps.begin(-20,-20,40,20,"balj");
		//ps.polygon(this->lattice_corners_,"linecolor=green");
		//ps.put(this->x_[0](0),this->x_[0](1),"0");

		Vector<double> x_loop(this->x_[0]);
		bool check_if_loop(false);
		for(unsigned int i(1);i<this->n_;i++){
			if(lattice_type_){
				if(i%2){ this->x_[i] = this->x_[i-1] + this->dir_nn_[0]; }
				else   { this->x_[i] = this->x_[i-1] - this->dir_nn_[1]; }
			} else {
				if(i%2){ this->x_[i] = this->x_[i-1] + this->dir_nn_[1]; }
				else   { this->x_[i] = this->x_[i-1] - this->dir_nn_[2]; }
			}
			if(reset_pos_in_lattice(this->x_[i])){ check_if_loop = true; }
			if(check_if_loop && my::are_equal(this->x_[i],x_loop)){
				check_if_loop = false;
				if(lattice_type_){ this->x_[i](1) += sqrt(3.0); }
				else { this->x_[i] += this->dir_nn_[0]-this->dir_nn_[1]; }
				reset_pos_in_lattice(this->x_[i]);
				x_loop = this->x_[i];
			}
			//this->x_[i] = this->x_[i].chop();
			//ps.put(this->x_[i](0),this->x_[i](1),my::tostring(i));
		}
		//ps.end(true,true,true);

		if(!this->obs_.size()){
			Vector<unsigned int> l(2);
			l(0) = 3;
			l(1) = 0;
			this->set_nn_links(l);
		}

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size()){
			if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
			else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
		}
	}
}

template<typename Type>
Honeycomb<Type>::~Honeycomb() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Honeycomb<Type>::set_obs(int nobs){
	if(nobs<0){ nobs = 1; }
	unsigned int nlinks;
	unsigned int nval;
	if(nobs>0){/*bond energy (valid for Honeycomb0pp*/
		nlinks = this->obs_[0].nlinks();
		nval = this->z_*this->spuc_/2;
		this->obs_.push_back(Observable("Bond energy",1,nval,nlinks));
		this->obs_[1].remove_links();
		for(unsigned int i(0);i<nlinks;i++){
			this->obs_[0](i,2) = this->get_site_in_ab(this->obs_[0](i,0))/2*3;
			this->obs_[0](i,2)+=(this->get_site_in_ab(this->obs_[0](i,1))-1)/2;
		}
	}
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Honeycomb<Type>::set_geometry(unsigned int const& n){
	L_ = sqrt(n/2.0);
	if(my::are_equal(L_,floor(L_))){
		lattice_type_ = 0;
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
		lattice_type_ = 1;
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
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry"<<std::endl;
	return Matrix<double>();
}

template<typename Type>
bool Honeycomb<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(lattice_type_){
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
