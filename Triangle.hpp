#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "System2DBis.hpp"

template<typename Type>
class Triangle: public System2DBis<Type>{
	public:
		/*!Constructor that organises the n=L^2 sites (L integer)*/
		Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Triangle()=0;

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
Triangle<Type>::Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2DBis<Type>(set_geometry(this->n_),ab,spuc,6,filename)
{
	if(this->status_==2){
		if(!this->obs_.size()){
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

			PSTricks ps("./","test");
			ps.begin(-20,-20,40,20,"balj");
			ps.put(this->x_[0](0),this->x_[0](1),"0");
			ps.polygon(this->lattice_corners_,"linecolor=green");

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
				ps.put(this->x_[i](0),this->x_[i](1),my::tostring(i));
			}
			ps.end(true,true,true);

			this->set_nn_links(Vector<unsigned int>(1,3));

			/*!sets the bond energy if it has not been set yet*/
			if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() == 1){
				this->J_.set(this->obs_[0].nlinks(),1);
			}
		}
	}
}

template<typename Type>
Triangle<Type>::~Triangle() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Triangle<Type>::set_obs(int nobs){
	(void)(nobs);
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Triangle<Type>::set_geometry(unsigned int const& n){
	L_ = sqrt(n/3.0);
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
		tmp(5,1) =-tmp(2,1);
		tmp(6,0) = tmp(0,0);
		tmp(6,1) = tmp(0,1);
		return tmp;
	}
	L_ = sqrt(n)/3.0;
	if(my::are_equal(L_,floor(L_))){
		lattice_type_ = 1;
		double a(sqrt(3.0)/2);
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
	return Matrix<double>();
}

template<typename Type>
bool Triangle<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(lattice_type_){
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
