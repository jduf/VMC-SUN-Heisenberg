#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "System2D.hpp"

template<typename Type>
class Square: public System2D<Type>{
	public:
		/*!Constructor that organises the n=p^2+q^2 sites on a (tilted-)square lattice*/
		Square(unsigned int const& spuc, unsigned int const& length, unsigned int const& tilt, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Square()=0;

	protected:
		void set_obs(int nobs);
		/*!Returns the neighbours of site i*/
		Vector<double> get_pos_in_lattice(unsigned int const& i) const;

	private:
		Matrix<double> set_ab(unsigned int const& spuc, unsigned int const& length, unsigned int const& tilt) const;
		Matrix<double> set_geometry(unsigned int const& n) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const;
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const;
		Vector<double> set_linear_jump() const;
};

/*{constructor*/
template<typename Type>
Square<Type>::Square(unsigned int const& spuc, unsigned int const& length, unsigned int const& tilt, std::string const& filename):
	System2D<Type>(set_geometry(this->n_),set_ab(spuc,length,tilt),set_linear_jump(),spuc,4,filename)
{
	if(this->status_==2){
		if(!this->obs_.size()){
			this->dir_nn_LxLy_.set(4,2);

			Vector<double> dir(2);
			dir(0) = 1.0;
			dir(1) = 0.0;
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(0,0) = dir(0);
			this->dir_nn_LxLy_(0,1) = dir(1);

			dir(0) = 0.0;
			dir(1) = 1.0;
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(1,0) = dir(0);
			this->dir_nn_LxLy_(1,1) = dir(1);

			dir(0) =-1.0;
			dir(1) = 0.0;
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(2,0) = dir(0);
			this->dir_nn_LxLy_(2,1) = dir(1);

			dir(0) = 0.0;
			dir(1) =-1.0;
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(3,0) = dir(0);
			this->dir_nn_LxLy_(3,1) = dir(1);

			this->set_nn_links(Vector<unsigned int>(1,2));
		}

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size()){
			if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
			else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
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
		this->obs_.push_back(Observable(this->n_,this->n_,50,5,false));
		for(unsigned int i(0);i<this->n_;i++){
			this->obs_[2](i,0) = 0;
			this->obs_[2](i,1) = i;
			this->obs_[2](i,2) = i;
		}
	}
}

template<typename Type>
Vector<double> Square<Type>::get_pos_in_lattice(unsigned int const& i) const {
	Vector<double> tmp(2);
	tmp(0) = i;
	tmp(1) = i/this->xloop_;
	this->set_pos_LxLy(tmp);
	return (this->LxLy_*tmp).chop();
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Square<Type>::set_ab(unsigned int const& spuc, unsigned int const& length, unsigned int const& tilt) const {
	if(!length){ return set_geometry(spuc); }
	if(!(spuc%length)){
		Matrix<double> tmp(2,2);
		tmp(0,0) = length;
		tmp(1,0) = 0;
		tmp(0,1) = tilt;
		tmp(1,1) = spuc/length;
		return tmp;
	} 
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown unit cell geometry"<<std::endl;
	return Matrix<double>();
}

template<typename Type>
Matrix<double> Square<Type>::set_geometry(unsigned int const& n) const {
	for(unsigned int p(0);p<=sqrt(n);p++){
		for(unsigned int q(0);q<p+1;q++){
			if(p*p+q*q==n){
				Matrix<double> tmp(2,2);
				tmp(0,0) = p;
				tmp(1,0) = q;
				tmp(0,1) = q;
				tmp(1,1) = (q?-double(p):p);
				return tmp;
			}
		}
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown lattice geometry, possible options are :"<<std::endl;
	for(unsigned int i(4);i<2*n;i++){
		for(unsigned int p(0);p<=sqrt(i);p++){
			for(unsigned int q(0);q<p+1;q++){
				if(p*p+q*q==i){
					std::cerr<<"(p="<<p<<", q="<<q<<") : n="<<i<<std::endl;
				}
			}
		}
	}
	return Matrix<double>();
}

template<typename Type>
Vector<double> Square<Type>::vector_towards(unsigned int const& i, unsigned int const& dir) const {
	(void)(i);
	Vector<double> tmp(2);
	tmp(0) = this->dir_nn_LxLy_(dir,0);
	tmp(1) = this->dir_nn_LxLy_(dir,1);
	return tmp;
}

template<typename Type>
void Square<Type>::try_neighbourg(Vector<double>& tn, unsigned int const& j) const {
	if(j%this->xloop_==0){
		tn(0) += this->dir_nn_LxLy_(1,0);
		tn(1) += this->dir_nn_LxLy_(1,1);
	}
	tn(0) += this->dir_nn_LxLy_(0,0);
	tn(1) += this->dir_nn_LxLy_(0,1);
}

template<typename Type>
Vector<double> Square<Type>::set_linear_jump() const {
	/*{!As the minimal unit cell contains 1 sites, the linear size is given by
	 * going from 0 to 1
	 *
	 * |  |
	 * 0--1--
	 *}*/
	Vector<double> tmp(2);
	tmp(0) = 1.0;
	tmp(1) = 0.0;
	return tmp;
}
/*}*/
#endif
