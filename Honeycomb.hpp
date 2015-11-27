#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "System2D.hpp"

template<typename Type>
class Honeycomb: public System2D<Type>{
	public:
		/*!Constructor that organises the n=2L^2 sites (L integer)*/
		Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Honeycomb()=0;

	protected:
		Matrix<double> dir_nn_;

		void set_observables(int nobs);
		/*!Returns the neighbours of site i*/
		Vector<double> get_pos_in_lattice(unsigned int const& i) const;

	private:
		Matrix<double> set_geometry(unsigned int const& n) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const;
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const;
		Vector<double> set_linear_jump() const;
};

/*{constructor*/
template<typename Type>
Honeycomb<Type>::Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(this->n_),ab,set_linear_jump(),spuc,3,filename),
	dir_nn_(3,2)
{
	if(this->status_==2){
		if(!this->obs_.size()){
			this->dir_nn_LxLy_.set(3,2);
			/*!as the linear_jump goes over two sites, xloop_ is twice bigger*/
			this->xloop_ *= 2;
			/*{!the directions are given for the sublattice with even site number
			 *  in the cartesian basis
			 * 
			 * (-1,sqrt(3))/2
			 *       \
			 *        x--(1,0)
			 *       / 
			 * (-1,-sqrt(3))/2
			 *
			 *  x = 0,2,4,...
			 *}*/
			Vector<double> dir(2);
			dir(0) = 1.0;
			dir(1) = 0.0;
			dir_nn_(0,0) = dir(0);
			dir_nn_(0,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(0,0) = dir(0);
			this->dir_nn_LxLy_(0,1) = dir(1);

			dir(0) =-0.5;
			dir(1) = sqrt(3.0)/2.0;
			dir_nn_(1,0) = dir(0);
			dir_nn_(1,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(1,0) = dir(0);
			this->dir_nn_LxLy_(1,1) = dir(1);

			dir(0) =-0.5;
			dir(1) =-sqrt(3.0)/2.0;
			dir_nn_(2,0) = dir(0);
			dir_nn_(2,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(2,0) = dir(0);
			this->dir_nn_LxLy_(2,1) = dir(1);

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
void Honeycomb<Type>::set_observables(int nobs){
	this->E_.set(50,5,false);

	if(nobs<0){ nobs = 1; }
	if(nobs>1){ /*the long range correlation*/
		this->obs_.push_back(Observable(this->n_,this->n_,50,5,false));
		for(unsigned int i(0);i<this->n_;i++){
			this->obs_[1](i,0) = 0;
			this->obs_[1](i,1) = i;
			this->obs_[1](i,2) = i;
		}
	}
}

template<typename Type>
Vector<double> Honeycomb<Type>::get_pos_in_lattice(unsigned int const& i) const {
	unsigned int j(i/2);
	Vector<double> tmp(this->linear_jump_*j);
	j = i/this->xloop_;
	tmp(0) += j*(dir_nn_(0,0)-dir_nn_(2,0));
	tmp(1) += j*(dir_nn_(0,1)-dir_nn_(2,1));
	j = i%this->xloop_;
	if(j%2){
		tmp(0) += dir_nn_(0,0);
		tmp(1) += dir_nn_(0,1);
	}
	this->set_pos_LxLy(tmp);
	return (this->LxLy_*tmp).chop();
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Honeycomb<Type>::set_geometry(unsigned int const& n) const {
	double L(sqrt(n/2));
	if(my::are_equal(L,floor(L))){
		Matrix<double> tmp(2,2);
		tmp(0,0) = L*1.5;
		tmp(1,0) = L*sqrt(3.0)/2.0;
		tmp(0,1) = L*1.5;
		tmp(1,1) =-L*sqrt(3.0)/2.0;
		return tmp;
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry"<<std::endl; 
	for(unsigned int i(4);i<2*n;i+=2){
		for(unsigned int p(0);p<=sqrt(i/2);p++){
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
Vector<double> Honeycomb<Type>::vector_towards(unsigned int const& i, unsigned int const& dir) const {
	/*{!the directions were defined as follows
	 * 1
	 *  \
	 *   i--0
	 *  / 
	 * 2 
	 *}*/
	Vector<double> tmp(2);
	if(i%2){//i=1,3,5,..
		/*{!therefore the for odd sites it gives
		 *     2
		 *     /
		 * 0--i
		 *     \
		 *      1 
		 *}*/
		tmp(0) = -this->dir_nn_LxLy_(dir,0);
		tmp(1) = -this->dir_nn_LxLy_(dir,1);
	} else {
		tmp(0) = this->dir_nn_LxLy_(dir,0);
		tmp(1) = this->dir_nn_LxLy_(dir,1);
	}
	return tmp;
}

template<typename Type>
void Honeycomb<Type>::try_neighbourg(Vector<double>& tn, unsigned int const& j) const {
	if(j%this->xloop_){
		if(j%2){//j=1,3,5,...
			tn(0) += this->dir_nn_LxLy_(0,0);
			tn(1) += this->dir_nn_LxLy_(0,1);
		} else {
			tn(0) -= this->dir_nn_LxLy_(1,0);
			tn(1) -= this->dir_nn_LxLy_(1,1);
		}
	} else {
		tn(0) += 2*this->dir_nn_LxLy_(0,0);
		tn(1) += 2*this->dir_nn_LxLy_(0,1);
	}
}

template<typename Type>
Vector<double> Honeycomb<Type>::set_linear_jump() const {
	/*{!As the minimal unit cell contains 2 sites, the linear size is given by
	 * going from 0 to 2
	 *
	 * 0--1
	 *     \
	 *      2
	 *}*/
	Vector<double> tmp(2);
	tmp(0) = 1.5;
	tmp(1) =-sqrt(3.0)/2.0;
	return tmp;
}
/*}*/
#endif
