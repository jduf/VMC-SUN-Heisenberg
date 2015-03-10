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

	protected:
		static Matrix<double> set_LxLy(unsigned int const& n){
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
};

template<typename Type>
Triangle<Type>::Triangle(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Triangle<Type>::set_LxLy(this->n_),ab,spuc,6,filename)
{
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : need to set the boundary condition"<<std::endl;
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : and need to check that they are correct"<<std::endl;
	if(this->status_==2){ this->compute_links(); }

	Vector<double> dir(2);
	dir(0) = 1.0;
	dir(1) = 0.0;
	dir = this->get_LxLy_pos(dir);
	this->dir_nn_LxLy_(0,0) = dir(0);
	this->dir_nn_LxLy_(0,1) = dir(1);

	dir(0) = 0.0;
	dir(1) = 1.0;
	dir = this->get_LxLy_pos(dir);
	this->dir_nn_LxLy_(1,0) = dir(0);
	this->dir_nn_LxLy_(1,1) = dir(1);

	dir(0) = -1.0;
	dir(1) = 1.0;
	dir = this->get_LxLy_pos(dir);
	this->dir_nn_LxLy_(2,0) = dir(0);
	this->dir_nn_LxLy_(2,1) = dir(1);

	dir(0) = -1.0;
	dir(1) = 0.0;
	dir = this->get_LxLy_pos(dir);
	this->dir_nn_LxLy_(3,0) = dir(0);
	this->dir_nn_LxLy_(3,1) = dir(1);

	dir(0) = 0.0;
	dir(1) = -1.0;
	dir = this->get_LxLy_pos(dir);
	this->dir_nn_LxLy_(4,0) = dir(0);
	this->dir_nn_LxLy_(4,1) = dir(1);

	dir(0) = 1.0;
	dir(1) = -1.0;
	dir = this->get_LxLy_pos(dir);
	this->dir_nn_LxLy_(5,0) = dir(0);
	this->dir_nn_LxLy_(5,1) = dir(1);

	std::cout<<"arg"<<std::endl;

	std::cout<<this->dir_nn_LxLy_<<std::endl;

}

template<typename Type>
Triangle<Type>::~Triangle(){}
#endif
