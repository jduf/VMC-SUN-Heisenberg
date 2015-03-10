#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "System2D.hpp"

template<typename Type>
class Square: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(4,filename), to construct a system with 4 links
		 * per sites */
		/*}*/
		Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Square()=0;

	private:
		static Matrix<double> set_LxLy(unsigned int const& n) {
			Matrix<double> tmp;
			for(unsigned int p(0);p<=sqrt(n);p++){
				for(unsigned int q(0);q<p+1;q++){
					if(p*p+q*q==n){ 
						tmp.set(2,2);
						tmp(0,0) = p;
						tmp(1,0) = q;
						tmp(0,1) = q;
						tmp(1,1) = (q?-double(p):p);
						return tmp;
					}
				}
			}
			return tmp;
		}
};

/*{constructor*/
template<typename Type>
Square<Type>::Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Square<Type>::set_LxLy(this->n_),ab,spuc,4,filename)
{
	this->status_=2;
	if(this->status_==2){ 

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
		dir(1) = 0.0;
		dir = this->get_LxLy_pos(dir);
		this->dir_nn_LxLy_(2,0) = dir(0);
		this->dir_nn_LxLy_(2,1) = dir(1);

		dir(0) = 0.0;
		dir(1) = -1.0;
		dir = this->get_LxLy_pos(dir);
		this->dir_nn_LxLy_(3,0) = dir(0);
		this->dir_nn_LxLy_(3,1) = dir(1);
		std::cout<<"arg"<<std::endl;

		std::cout<<this->dir_nn_LxLy_<<std::endl;


		this->compute_links(); 
	}
}

template<typename Type>
Square<Type>::~Square(){}
/*}*/

/*{protected methods*/
/*}*/
#endif
