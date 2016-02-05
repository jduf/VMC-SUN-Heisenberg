#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "System2DBis.hpp"

template<typename Type>
class Triangle: public System2DBis<Type>{
	public:
		/*!Constructor*/
		Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Triangle()=0;

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
Triangle<Type>::Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2DBis<Type>(set_geometry((!this->obs_.size() || !this->obs_[0].nlinks())?this->n_:0),ab,spuc,6,filename)
{}

template<typename Type>
Triangle<Type>::~Triangle() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Triangle<Type>::init_lattice(){
	if(!this->obs_.size() || !this->obs_[0].nlinks()){
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

		this->x_[0](0)+= 0.01;
		this->x_[0](1)+= 0.01;

		Vector<double> x_loop(this->x_[0]);
		bool check_if_loop(false);
		for(unsigned int i(1);i<this->n_;i++){
			this->x_[i] = this->x_[i-1] + this->dir_nn_[0];
			if(reset_pos_in_lattice(this->x_[i])){ check_if_loop = true; }
			if(check_if_loop && my::are_equal(this->x_[i],x_loop)){
				check_if_loop = false;
				this->x_[i] += this->dir_nn_[1];
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

		this->set_nn_links(Vector<unsigned int>(1,3));

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() == 1){
			this->J_.set(this->obs_[0].nlinks(),1);
		}
	}
}

template<typename Type>
void Triangle<Type>::set_obs(int nobs){
	(void)(nobs);
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Triangle<Type>::set_geometry(unsigned int const& n){
	if(n){
		L_ = sqrt(n/3.0);
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
		L_ = sqrt(n)/3.0;
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2.0);
			Matrix<double> tmp(7,2);
			tmp(0,0) = 0.0;
			tmp(0,1) =-2.0*a*L_;
			tmp(1,0) = 1.5*L_;
			tmp(1,1) =-a*L_;
			tmp(2,0) = tmp(1,0);
			tmp(2,1) =-tmp(1,1);
			tmp(3,0) = tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(2,0);
			tmp(5,1) =-tmp(2,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}
		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		for(unsigned int l(2);l<20;l++){ std::cerr<<"n="<<3*l*l<<" or "<<(3*l)*(3*l)<<std::endl; }
		std::cerr<<"n=3*l*l or (3*l)^2"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Triangle<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(this->ref_(3)){
			double t(tan(M_PI/6.0)*x(0)/x(1));
			if(x(0)>0){
				if(std::abs(t)>1){ x+=this->dir_nn_[3]*L_*3.0; }
				else {
					if(t>0){       x+=this->dir_nn_[4]*L_*3.0; }
					else   {       x+=this->dir_nn_[2]*L_*3.0; }
				}
			} else {
				if(std::abs(t)>1){ x+=this->dir_nn_[0]*L_*3.0; }
				else {
					if(t>0){       x+=this->dir_nn_[1]*L_*3.0; }
					else   {       x+=this->dir_nn_[5]*L_*3.0; }
				}
			}
		} else {
			double t(tan(M_PI/3.0)*x(0)/x(1));
			if(x(1)>0){
				if(std::abs(t)<1){ x+=(this->dir_nn_[4]+this->dir_nn_[5])*L_; }
				else {
					if(t>0){       x+=(this->dir_nn_[4]+this->dir_nn_[3])*L_; }
					else   {       x+=(this->dir_nn_[5]+this->dir_nn_[0])*L_; }
				}
			} else {
				if(std::abs(t)<1){ x+=(this->dir_nn_[1]+this->dir_nn_[2])*L_; }
				else {
					if(t>0){       x+=(this->dir_nn_[1]+this->dir_nn_[0])*L_; }
					else   {       x+=(this->dir_nn_[2]+this->dir_nn_[3])*L_; }
				}
			}
		}
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Triangle<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const {
	(void)(i);
	return this->dir_nn_[d];
}
/*}*/
#endif
