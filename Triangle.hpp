#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "System2D.hpp"

template<typename Type>
class Triangle: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(6,filename), to construct a system with 6 links
		 * per sites */
		/*}*/
		Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Triangle()=0;

		Vector<double> compute_J(Vector<double> const& Jp);

	protected:
		Vector<double> get_pos_in_lattice(unsigned int const& i) const;

	private:
		Matrix<double> set_LxLy(unsigned int const& n) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const;
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const;
};

/*{constructor*/
template<typename Type>
Triangle<Type>::Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Triangle<Type>::set_LxLy(this->n_),ab,spuc,6,filename)
{
	if(this->status_==2){ 
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

		dir(0) = -1.0;
		dir(1) = 1.0;
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(2,0) = dir(0);
		this->dir_nn_LxLy_(2,1) = dir(1);

		dir(0) = -1.0;
		dir(1) = 0.0;
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(3,0) = dir(0);
		this->dir_nn_LxLy_(3,1) = dir(1);

		dir(0) = 0.0;
		dir(1) = -1.0;
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(4,0) = dir(0);
		this->dir_nn_LxLy_(4,1) = dir(1);

		dir(0) = 1.0;
		dir(1) = -1.0;
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(5,0) = dir(0);
		this->dir_nn_LxLy_(5,1) = dir(1);

		//this->compute_nn_links(); 
	}
}

template<typename Type>
Triangle<Type>::~Triangle() = default;
/*}*/

/*{protected methods*/
template<typename Type>
Vector<double> Triangle<Type>::get_pos_in_lattice(unsigned int const& i) const {
	Vector<double> tmp(2);
	tmp(0) = i;
	tmp(1) = i/this->xloop_;
	return tmp;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Triangle<Type>::set_LxLy(unsigned int const& n) const {
	Matrix<double> tmp;
	if(n==27){
		tmp.set(2,2);
		tmp(0,0) = 6;
		tmp(1,0) = -3;
		tmp(0,1) = 3;
		tmp(1,1) = 3;
	}
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
/*}*/
#endif

