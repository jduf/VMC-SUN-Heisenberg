#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "System2D.hpp"

/*{Description*/
/*!the allowed clusters are that n=2L^2, L integer.
 * as the minimal unit cell contains two sites, the lattice basis vectors are
 * given by the ones that are jumping from one horizontal pair to the
 * neighbouring ones. in the orthogonal basis where the bound lenght is one,
 * this lattice basis is given by ((1/3,1/3),(-1/2,1/2))
 *
 * the cluster basis vectors (Lx,Ly) are colinear with the lattice ones, so one
 * can express the cluster basis as ((L,0),(0,L))
 *
 * the unit cell basis has to be expressed with the lattice basis vectors
 */
/*}*/
template<typename Type>
class Honeycomb: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises*/
		/*}*/
		Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Honeycomb()=0;

	protected:
		void set_observables(int nobs);
		/*!Returns the neighbours of site i*/
		Vector<double> get_pos_in_lattice(unsigned int const& i) const;

	private:
		Matrix<double> dir_nn_;

		Matrix<double> set_geometry(unsigned int const& n) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const;
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const;
};

/*{constructor*/
template<typename Type>
Honeycomb<Type>::Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(this->n_),ab,spuc,3,filename),
	dir_nn_(this->z_,2)
{
	if(this->status_==2){
		this->xloop_ *= 2;
		/*!the directions are given for the sublattice with even site number
		 * in the cartesian basis
		 * 
		 * (-1,sqrt(3))/2
		 *       \
		 *        x--(1,0)
		 *       / 
		 * (-1,-sqrt(3))/2
		 *
		 *  x = 0,2,4,...
		 * */
		Vector<double> dir(2);
		dir(0) = 1.0;
		dir(1) = 0.0;
		dir_nn_(0,0) = dir(0);
		dir_nn_(0,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(0,0) = dir(0);
		this->dir_nn_LxLy_(0,1) = dir(1);

		dir(0) = -0.5;
		dir(1) = sqrt(3.0)/2.0;
		dir_nn_(1,0) = dir(0);
		dir_nn_(1,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(1,0) = dir(0);
		this->dir_nn_LxLy_(1,1) = dir(1);

		dir(0) = -0.5;
		dir(1) = -sqrt(3.0)/2.0;
		dir_nn_(2,0) = dir(0);
		dir_nn_(2,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(2,0) = dir(0);
		this->dir_nn_LxLy_(2,1) = dir(1);

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
Honeycomb<Type>::~Honeycomb() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Honeycomb<Type>::set_observables(int nobs){
	this->E_.set(50,5,false);
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
	Vector<double> tmp(2);
	unsigned int j(i/this->xloop_);
	tmp(0) = j*(dir_nn_(0,0)-dir_nn_(2,0));
	tmp(1) = j*(dir_nn_(0,1)-dir_nn_(2,1));
	j = i-j*this->xloop_;
	tmp(0) += j/2*(dir_nn_(0,0)-dir_nn_(1,0));
	tmp(1) += j/2*(dir_nn_(0,1)-dir_nn_(1,1));
	if(j%2==1){
		tmp(0) += dir_nn_(0,0);
		tmp(1) += dir_nn_(0,1);
	}
	return tmp;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Honeycomb<Type>::set_geometry(unsigned int const& n) const {
	Matrix<double> tmp;
	double L(sqrt(n/2));
	if(my::are_equal(L,floor(L))){
		tmp.set(2,2);
		tmp(0,0) = L*1.5;
		tmp(1,0) = L*sqrt(3.0)/2.0;
		tmp(0,1) = L*1.5;
		tmp(1,1) =-L*sqrt(3.0)/2.0;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : unkown geometry"<<std::endl; }
	return tmp;
}

template<typename Type>
Vector<double> Honeycomb<Type>::vector_towards(unsigned int const& i, unsigned int const& dir) const {
	Vector<double> tmp(2);
	if(i%2){//i=1,3,5,..
		tmp(0) = -this->dir_nn_LxLy_(dir,0);
		tmp(1) = -this->dir_nn_LxLy_(dir,1);
	} else {
		tmp(0) = this->dir_nn_LxLy_(dir,0);
		tmp(1) = this->dir_nn_LxLy_(dir,1);
	}
	return tmp;
}

template<typename Type>
void Honeycomb<Type>::try_neighbourg(Vector<double>& tn, unsigned int const& i) const {
	if(i%this->xloop_==0){
		tn(0) += 2*this->dir_nn_LxLy_(0,0);
		tn(1) += 2*this->dir_nn_LxLy_(0,1);
	} else {
		if(i%2){//i=1,3,5,...
			tn(0) += this->dir_nn_LxLy_(0,0);
			tn(1) += this->dir_nn_LxLy_(0,1);
		} else {
			tn(0) -= this->dir_nn_LxLy_(1,0);
			tn(1) -= this->dir_nn_LxLy_(1,1);
		}
	}
}
/*}*/
#endif
