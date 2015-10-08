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
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(3,filename), to construct a system with 3 links
		 * per sites */
		/*}*/
		Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Honeycomb()=0;

	protected:
		Vector<double> get_pos_in_lattice(unsigned int const& i) const;

	private:
		Matrix<double> dir_nn_;

		Matrix<double> set_LxLy(unsigned int const& n) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const;
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const;
};

/*{constructor*/
template<typename Type>
Honeycomb<Type>::Honeycomb(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Honeycomb<Type>::set_LxLy(this->n_),ab,spuc,3,filename),
	dir_nn_(this->z_,2)
{
	if(this->status_==2){
		this->xloop_ *= 2;
		/*the directions are given for the sublattice with odd site number*/
		double h(1./3.);
		Vector<double> dir(2);
		dir(0) = h;
		dir(1) = h;
		dir_nn_(0,0) = dir(0);
		dir_nn_(0,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(0,0) = dir(0);
		this->dir_nn_LxLy_(0,1) = dir(1);

		dir(0) = h-1.0;
		dir(1) = h;
		dir_nn_(1,0) = dir(0);
		dir_nn_(1,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(1,0) = dir(0);
		this->dir_nn_LxLy_(1,1) = dir(1);

		dir(0) = h;
		dir(1) = h-1.0;
		dir_nn_(2,0) = dir(0);
		dir_nn_(2,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(2,0) = dir(0);
		this->dir_nn_LxLy_(2,1) = dir(1);

		std::cerr<<__PRETTY_FUNCTION__<<" : new def of set_nn_links will be problematic"<<std::endl;
		this->set_nn_links(Vector<unsigned int>(1,3));
	}
}

template<typename Type>
Honeycomb<Type>::~Honeycomb() = default;
/*}*/

/*{protected methods*/
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
Matrix<double> Honeycomb<Type>::set_LxLy(unsigned int const& n) const {
	Matrix<double> tmp;
	double L(sqrt(n/2));
	if(my::are_equal(L,floor(L))){
		tmp.set(2,2);
		tmp(0,0) = L;
		tmp(1,0) = 0.0;
		tmp(0,1) = 0.0;
		tmp(1,1) = L;
	}
	return tmp;
}

template<typename Type>
Vector<double> Honeycomb<Type>::vector_towards(unsigned int const& i, unsigned int const& dir) const {
	Vector<double> tmp(2);
	if(i%2==0){
		tmp(0) = this->dir_nn_LxLy_(dir,0);
		tmp(1) = this->dir_nn_LxLy_(dir,1);
	} else {
		tmp(0) = -this->dir_nn_LxLy_((dir+1)%3,0);
		tmp(1) = -this->dir_nn_LxLy_((dir+1)%3,1);
	}
	return tmp;
}

template<typename Type>
void Honeycomb<Type>::try_neighbourg(Vector<double>& tn, unsigned int const& j) const {
	if(j%2==0){
		tn(0) -= this->dir_nn_LxLy_(1,0);
		tn(1) -= this->dir_nn_LxLy_(1,1);
	} else {
		tn(0) += this->dir_nn_LxLy_(0,0);
		tn(1) += this->dir_nn_LxLy_(0,1);
	}
	if(j%this->xloop_==0){
		tn(0) += this->dir_nn_LxLy_(0,0);
		tn(1) += this->dir_nn_LxLy_(0,1);
		tn(0) -= this->dir_nn_LxLy_(2,0);
		tn(1) -= this->dir_nn_LxLy_(2,1);
	}
}
/*}*/

//template<typename Type>
//Matrix<int> Honeycomb<Type>::get_neighbourg(unsigned int i) const {
//Matrix<int> nb(this->z_,2,1);
//switch(this->spuc_){
//case 6:
//{
//switch(i%6){
//case 0:
//{
///*+x+y neighbour*/
//nb(0,0) = i+1;
///*-x neighbour*/
//if(i%(this->spuc_*this->Lx_)){ nb(1,0) = i-3; }
//else {
//nb(1,0) = i-3+this->spuc_*this->Lx_;
//nb(1,1) = this->bc_;
//}
///*+x-y neighbour*/
//nb(2,0) = i+5;
//}break;
//case 1:
//{
///*+x neighbour*/
//nb(0,0) = i+1;
///*-x+y neighbour*/
//if(i<this->n_-this->spuc_*this->Lx_){
//if((i-1)%(this->spuc_*this->Lx_)){ nb(1,0) = i-3+this->spuc_*this->Lx_ ; }
//else {
//nb(1,0) = i-3+2*this->spuc_*this->Lx_;
//nb(1,1) = this->bc_;
//}
//} else {
//if(i-1!=this->n_-this->spuc_*this->Lx_){
//nb(1,0) = i-3-this->spuc_*this->Lx_*(this->Ly_-1);
//nb(1,1) = this->bc_;
//} else {
//nb(1,0) = -2+this->spuc_*this->Lx_;
//nb(1,1) = this->bc_*this->bc_;
//}
//}
///*-x-y neighbour*/
//nb(2,0) = i-1;
//}break;
//case 2:
//{
///*+x-y neighbour*/
//nb(0,0) = i+1;
///*+x+y neighbour*/
//if(i<this->n_-this->spuc_*this->Lx_){ nb(1,0) = i+3+this->spuc_*this->Lx_; }
//else {
//nb(1,0) = i+3-this->spuc_*this->Lx_*(this->Ly_-1);
//nb(1,1) = this->bc_;
//}
///*-x neighbour*/
//nb(2,0) = i-1;
//}break;
//case 3:
//{
///*-x-y neighbour*/
//nb(0,0) = i+1;
///*+x neighbour*/
//if((i+3)%(this->spuc_*this->Lx_)){ nb(1,0) = i+3; }
//else {
//nb(1,0) = i+3-this->spuc_*this->Lx_;
//nb(1,1) = this->bc_;
//}
///*-x+y neighbour*/
//nb(2,0) = i-1;
//}break;
//case 4:
//{
///*-x neighbour*/
//nb(0,0) = i+1;
///*+x-y neighbour*/
//if(i>this->spuc_*this->Lx_){
//if((i+2)%(this->spuc_*this->Lx_)){ nb(1,0) = i+3-this->spuc_*this->Lx_; }
//else {
//nb(1,0) = i+3-2*this->spuc_*this->Lx_;
//nb(1,1) = this->bc_;
//}
//} else {
//if(i+2!=this->spuc_*this->Lx_){
//nb(1,0) = i+3+this->spuc_*this->Lx_*(this->Ly_-1);
//nb(1,1) = this->bc_;
//} else {
//nb(1,0) = this->n_+1-this->spuc_*this->Lx_;
//nb(1,1) = this->bc_*this->bc_;
//}
//}
///*+x+y neighbour*/
//nb(2,0) = i-1;
//}break;
//case 5:
//{
///*-x+y neighbour*/
//nb(0,0) = i-5;
///*-x-y neighbour*/
//if(i>this->spuc_*this->Lx_){ nb(1,0) = i-3-this->spuc_*this->Lx_; }
//else {
//nb(1,0) = i-3+this->spuc_*this->Lx_*(this->Ly_-1);
//nb(1,1) = this->bc_;
//}
///*+x neighbour*/
//nb(2,0) = i-1;
//}break;
//}
//}break;
//}
//return nb;
//}
#endif
