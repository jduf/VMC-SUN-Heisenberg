#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "System2D.hpp"

template<typename Type>
class Triangle: public System2D<Type>{
	public:
		/*!Constructor that organises the n=L^2 sites (L integer)*/
		Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Triangle()=0;

	protected:
		Matrix<double> dir_nn_;

		void set_obs(int nobs);
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
Triangle<Type>::Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(this->n_),ab,set_linear_jump(),spuc,6,filename),
	dir_nn_(6,2)
{
	if(this->status_==2){
		if(!this->obs_.size()){
			this->dir_nn_LxLy_.set(6,2);

			Vector<double> dir(2);
			dir(0) = 1.0;
			dir(1) = 0.0;
			dir_nn_(0,0) = dir(0);
			dir_nn_(0,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(0,0) = dir(0);
			this->dir_nn_LxLy_(0,1) = dir(1);

			dir(0) = 0.5;
			dir(1) = sqrt(3.0)/2.0;
			dir_nn_(1,0) = dir(0);
			dir_nn_(1,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(1,0) = dir(0);
			this->dir_nn_LxLy_(1,1) = dir(1);

			dir(0) =-0.5;
			dir(1) = sqrt(3.0)/2.0;
			dir_nn_(2,0) = dir(0);
			dir_nn_(2,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(2,0) = dir(0);
			this->dir_nn_LxLy_(2,1) = dir(1);

			dir(0) =-1.0;
			dir(1) = 0.0;
			dir_nn_(3,0) = dir(0);
			dir_nn_(3,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(3,0) = dir(0);
			this->dir_nn_LxLy_(3,1) = dir(1);

			dir(0) =-0.5;
			dir(1) =-sqrt(3.0)/2.0;
			dir_nn_(4,0) = dir(0);
			dir_nn_(4,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(4,0) = dir(0);
			this->dir_nn_LxLy_(4,1) = dir(1);

			dir(0) = 0.5;
			dir(1) =-sqrt(3.0)/2.0;
			dir_nn_(5,0) = dir(0);
			dir_nn_(5,1) = dir(1);
			this->set_pos_LxLy(dir);
			this->dir_nn_LxLy_(5,0) = dir(0);
			this->dir_nn_LxLy_(5,1) = dir(1);

			this->set_nn_links(Vector<unsigned int>(1,3));
		}

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() == 1){
			this->J_.set(this->obs_[0].nlinks(),1);
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

template<typename Type>
Vector<double> Triangle<Type>::get_pos_in_lattice(unsigned int const& i) const {
	Vector<double> tmp(this->linear_jump_*i);
	unsigned int j(i/this->xloop_);
	tmp(0) += j*dir_nn_(1,0);
	tmp(1) += j*dir_nn_(1,1);
	this->set_pos_LxLy(tmp);
	return this->LxLy_*tmp;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Triangle<Type>::set_geometry(unsigned int const& n) const {
	Matrix<double> tmp;
	double L(sqrt(n));
	if(my::are_equal(L,floor(L)) && int(L)%2==0){
		tmp.set(2,2);
		tmp(0,0) = L;
		tmp(1,0) = 0;
		tmp(0,1) = L/2;
		tmp(1,1) = L*sqrt(3.0)/2.0;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry"<<std::endl; }
	return tmp;
}

template<typename Type>
Vector<double> Triangle<Type>::vector_towards(unsigned int const& i, unsigned int const& dir) const {
	(void)(i);
	Vector<double> tmp(2);
	tmp(0) = this->dir_nn_LxLy_(dir,0);
	tmp(1) = this->dir_nn_LxLy_(dir,1);
	return tmp;
}

template<typename Type>
void Triangle<Type>::try_neighbourg(Vector<double>& tn, unsigned int const& j) const {
	if(j%this->xloop_==0){
		tn(0) += this->dir_nn_LxLy_(1,0);
		tn(1) += this->dir_nn_LxLy_(1,1);
	}
	tn(0) += this->dir_nn_LxLy_(0,0);
	tn(1) += this->dir_nn_LxLy_(0,1);
}

template<typename Type>
Vector<double> Triangle<Type>::set_linear_jump() const {
	Vector<double> tmp(2);
	tmp(0) = 1.0;
	tmp(1) = 0.0;
	return tmp;
}
/*}*/
#endif
