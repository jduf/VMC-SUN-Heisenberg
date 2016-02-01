#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "System2DBis.hpp"

template<typename Type>
class Square: public System2DBis<Type>{
	public:
		/*!Constructor*/
		Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Square()=0;

	protected:
		void set_obs(int nobs);

	private:
		unsigned int p_;
		unsigned int q_;

		Matrix<double> set_geometry(unsigned int const& n);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const;
};

/*{constructor*/
template<typename Type>
Square<Type>::Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2DBis<Type>(set_geometry((!this->obs_.size() || !this->obs_[0].nlinks())?this->n_:0),ab,spuc,4,filename)
{
	if(this->status_==2){
		if(!this->obs_.size() || !this->obs_[0].nlinks()){
			/*{!the directions are given in the cartesian basis
			 *
			 *       (1,0)
			 *         |
			 * (-1,0)--x--(1,0)
			 *         | 
			 *      (-1,0)
			 *}*/
			this->dir_nn_[0](0) = 1.0;
			this->dir_nn_[0](1) = 0.0;

			this->dir_nn_[1](0) = 0.0;
			this->dir_nn_[1](1) = 1.0;

			this->dir_nn_[2](0) =-1.0;
			this->dir_nn_[2](1) = 0.0;

			this->dir_nn_[3](0) = 0.0;
			this->dir_nn_[3](1) =-1.0;

			this->x_[0](0) = 0.02;
			this->x_[0](1) = 0.01;

			PSTricks ps("./","test");
			ps.begin(-20,-20,40,20,"balj");
			ps.polygon(this->lattice_corners_,"linecolor=green");
			ps.put(this->x_[0](0),this->x_[0](1),"0");
			
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

			for(unsigned int i(0);i<4;i++){
				this->boundary_[i].set(2);
				this->boundary_[i](0) = this->lattice_corners_(i,0);
				this->boundary_[i](1) = this->lattice_corners_(i,1);
			}

			this->set_nn_links(Vector<unsigned int>(1,2));

			/*!sets the bond energy if it has not been set yet*/
			if(this->obs_[0].nlinks() != this->J_.size()){
				if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
				else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
			}
		}
	}
}

template<typename Type>
Square<Type>::~Square() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Square<Type>::set_obs(int nobs){
	if(nobs<0){ nobs = 1; }
	if(nobs>1){ /*the long range correlation*/
		/*bond energy missing*/
		this->obs_.push_back(Observable("Long range correlations",2,this->n_,this->n_));
		for(unsigned int i(0);i<this->n_;i++){
			this->obs_[2](i,0) = 0;
			this->obs_[2](i,1) = i;
			this->obs_[2](i,2) = i;
		}
	}
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Square<Type>::set_geometry(unsigned int const& n){
	if(n){
		for(unsigned int p(0);p<=sqrt(n);p++){
			for(unsigned int q(0);q<p+1;q++){
				if(p*p+q*q==n){
					p_ = sqrt(n);
					if(my::are_equal(sqrt(n),p_) && my::get_yn("Two possible clusters, use the tilded one ?")){
						p_ = p;
						q_ = q;
					} else { q_ = 0; }
					Matrix<double> tmp(5,2);
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
			}
		}
		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		for(unsigned int p(0);p<20;p++){
			for(unsigned int q(p);q<p+1;q++){
				std::cerr<<"n="<<p*p+q*q<<std::endl;
			}
		}
		std::cerr<<"n=p*p+q*q (p>=q)"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Square<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		double a(this->lattice_corners_(0,0)*x(0)+this->lattice_corners_(0,1)*x(1));
		double b(this->lattice_corners_(1,0)*x(0)+this->lattice_corners_(1,1)*x(1));
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
Vector<double> Square<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const {
	(void)(i);
	return this->dir_nn_[d];
}
/*}*/
#endif
