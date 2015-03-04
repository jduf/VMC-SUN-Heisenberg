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
		Square(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Square()=0;

	protected:
		unsigned int xloop_;
		Matrix<double> ab_basis_;

		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int i) const;
		/*!Returns the position of the site i in the basis (a,b)*/
		Vector<double> get_ab_pos(Vector<double> const& x) const;
		/*!Reset x so that it belongs to the square (a,b)*/
		bool set_in_zone(Vector<double>& x) const;
		
	private:
		Matrix<double> inv_basis_;
		Matrix<double> dir_ab_;

		/*!Set the neighbour of site i in direction dir in nb*/
		void find_neighbourg(unsigned int i, unsigned int dir, Matrix<int>& nb) const;
};

/*{constructor*/
template<typename Type>
Square<Type>::Square(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Lx,Ly,spuc,4,filename),
	xloop_(0),
	ab_basis_(2,2),
	inv_basis_(2,2),
	dir_ab_(4,2)
{
	this->status_=2;
	if(this->status_==2){ 
		for(unsigned int p(0);p<=sqrt(this->n_);p++){
			for(unsigned int q(0);q<p+1;q++){
				if(p*p+q*q==this->n_){ 
					ab_basis_(0,0) = p;
					ab_basis_(1,0) = q;
					ab_basis_(0,1) = q;
					ab_basis_(1,1) = (q?-double(p):p);
					inv_basis_(0,0) = (q?-double(p):p);
					inv_basis_(1,0) = -double(q);
					inv_basis_(0,1) = -double(q);
					inv_basis_(1,1) = p;
				}
			}
		}
		inv_basis_ /= ab_basis_(0,0)*ab_basis_(1,1)-ab_basis_(1,0)*ab_basis_(0,1);

		Vector<double> x(2);
		unsigned int j(0);
		do{
			x(1) = 0;
			x(0) = ++j;
			x = get_ab_pos(x);
		} while( !are_equal(x(0),0) || !are_equal(x(1),0)  );
		xloop_ = j;

		Vector<double> dir(2);
		dir(0) = 1.0;
		dir(1) = 0.0;
		dir = get_ab_pos(dir);
		dir_ab_(0,0) = dir(0);
		dir_ab_(0,1) = dir(1);

		dir(0) = 0.0;
		dir(1) = 1.0;
		dir = get_ab_pos(dir);
		dir_ab_(1,0) = dir(0);
		dir_ab_(1,1) = dir(1);

		dir(0) = -1.0;
		dir(1) = 0.0;
		dir = get_ab_pos(dir);
		dir_ab_(2,0) = dir(0);
		dir_ab_(2,1) = dir(1);

		dir(0) = 0.0;
		dir(1) = -1.0;
		dir = get_ab_pos(dir);
		dir_ab_(3,0) = dir(0);
		dir_ab_(3,1) = dir(1);

		this->compute_links(); 
	}
}

template<typename Type>
Square<Type>::~Square(){}
/*}*/

/*{protected methods*/
template<typename Type>
Matrix<int> Square<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	find_neighbourg(i,0,nb);
	find_neighbourg(i,1,nb);
	find_neighbourg(i,2,nb);
	find_neighbourg(i,3,nb);
	return nb;
}

template<typename Type>
Vector<double> Square<Type>::get_ab_pos(Vector<double> const& x) const {
	Vector<double> tmp(inv_basis_*x);
	double ip;
	tmp(0) = std::modf(tmp(0),&ip);
	tmp(1) = std::modf(tmp(1),&ip);
	if( are_equal(tmp(0),1) ) { tmp(0) = 0; }
	if( are_equal(tmp(1),1) ) { tmp(1) = 0; }
	return tmp;
}

template<typename Type>
bool Square<Type>::set_in_zone(Vector<double>& x) const {
	bool in_zone(false);
	double ip;
	x(0) = std::modf(x(0),&ip);
	if( x(0)<0 ){ x(0) += 1.0; in_zone = !in_zone; }
	if( are_equal(x(0),1) ) { x(0) = 0; in_zone = !in_zone; }
	if( ip>0 ) { in_zone = !in_zone;}

	x(1) = std::modf(x(1),&ip);
	if( x(1)<0 ){ x(1) += 1.0;  in_zone = !in_zone; }
	if( are_equal(x(1),1) ) { x(1) = 0; in_zone = !in_zone; }
	if( ip>0 ) { in_zone = !in_zone ;}
	return in_zone;
}
/*}*/

/*{private methods*/
template<typename Type>
void Square<Type>::find_neighbourg(unsigned int i, unsigned int dir, Matrix<int>& nb) const{
	Vector<double> nn_ab(2,0);/*nearest neighbour in the ab basis*/
	Vector<double> tn_ab(2,0);/*trial neighbour in the ab basis*/
	Vector<double> tn_s(2,0); /*trial neighbour in the square basis*/
	unsigned int j(0);

	nn_ab(0) = i;
	nn_ab(1) = i/xloop_;

	nn_ab = get_ab_pos(nn_ab);
	set_in_zone(nn_ab);
	nn_ab(0) += dir_ab_(dir,0);
	nn_ab(1) += dir_ab_(dir,1);
	if(set_in_zone(nn_ab)){ nb(dir,1) = this->bc_; }
	do{
		tn_s(0) = j;
		tn_ab=get_ab_pos(tn_s);
		set_in_zone(tn_ab);
		if((j+1)%xloop_==0){ tn_s(1)+=1; }
		j++;
	} while ( !are_equal(tn_ab,nn_ab) && j<this->n_+2 );
	nb(dir,0) = j-1;
}
/*}*/
#endif
