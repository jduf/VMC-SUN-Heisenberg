#ifndef DEF_KAGOME
#define DEF_KAGOME

#include "System2D.hpp"

template<typename Type>
class Kagome: public System2D<Type>{
	public:
		/*!Constructor that organises the n=3L^2 sites (L integer)*/
		Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Kagome()=0;

	protected:
		Matrix<double> dir_nn_;

		void set_observables(int nobs);
		/*!Returns the neighbours of site i*/
		Vector<double> get_pos_in_lattice(unsigned int const& i) const;

	private:
		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& spuc) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const;
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const;
		Vector<double> set_linear_jump() const;
};

/*{constructor*/
template<typename Type>
Kagome<Type>::Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(this->n_,spuc),ab,set_linear_jump(),spuc,4,filename),
	dir_nn_(3,2)
{
	if(this->status_==2){
		this->dir_nn_LxLy_.set(3,2);
		/*!as the linear_jump goes over tree sites, xloop_ is tree times bigger*/
		this->xloop_ *= 3;
		/*{!the directions are given for the sublattice with even site number
		 *  in the cartesian basis
		 * 
		 * f--g
		 *  \/  
		 *  d   
		 *  /\
		 * a--b
		 *
		 * a-b = (1,0)
		 * b-d = (-1,sqrt(3))/2
		 * a-d = (1,sqrt(3))/2
		 *}*/
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

		dir(0) = 0.5;
		dir(1) = sqrt(3.0)/2.0;
		dir_nn_(2,0) = dir(0);
		dir_nn_(2,1) = dir(1);
		this->set_pos_LxLy(dir);
		this->dir_nn_LxLy_(2,0) = dir(0);
		this->dir_nn_LxLy_(2,1) = dir(1);

		Vector<unsigned int> l(3,2);
		this->set_nn_links(l);

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size()){
			if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
			else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
		}
	}
}

template<typename Type>
Kagome<Type>::~Kagome() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Kagome<Type>::set_observables(int nobs){
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
Vector<double> Kagome<Type>::get_pos_in_lattice(unsigned int const& i) const {
	unsigned int j(i/3);
	Vector<double> tmp(this->linear_jump_*j);
	j = i/this->xloop_;
	tmp(0) += 2*j*dir_nn_(1,0);
	tmp(1) += 2*j*dir_nn_(1,1);
	j = i%this->xloop_;
	switch(j%3){
		case 1:
			{
				tmp(0) += dir_nn_(0,0);
				tmp(1) += dir_nn_(0,1);
			}break;
		case 2:
			{
				tmp(0) += dir_nn_(2,0);
				tmp(1) += dir_nn_(2,1);
			}break;
	}
	this->set_pos_LxLy(tmp);
	return this->LxLy_*tmp;
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Kagome<Type>::set_geometry(unsigned int const& n, unsigned int const& spuc) const {
	Matrix<double> tmp;
	if(n%3){ std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry"<<std::endl; }
	else {
		if(spuc==9){
			double L(sqrt(n/spuc));
			if(my::are_equal(L,floor(L))){
				tmp.set(2,2);
				tmp(0,0) = 9;
				tmp(1,0) = -3*sqrt(3.0);
				tmp(0,1) = 0;
				tmp(1,1) = 6*sqrt(3.0);
			}
			if(n!=81){ std::cerr<<__PRETTY_FUNCTION__<<" : can make it work but need some modification here"<<std::endl; }
		} else {
			double L(sqrt(n/3));
			if(my::are_equal(L,floor(L))){
				tmp.set(2,2);
				tmp(0,0) = L*2;
				tmp(1,0) = 0;
				tmp(0,1) = L;
				tmp(1,1) = L*sqrt(3.0);
			}
		}
	}
	return tmp;
}

template<typename Type>
Vector<double> Kagome<Type>::vector_towards(unsigned int const& i, unsigned int const& dir) const {
	/*{The directions were defined as follows
	 * f--g
	 *  \/ 
	 *  d  
	 *  /\
	 * a--b
	 *
	 * a-b = 0
	 * b-d = 1
	 * a-d = 2 
	 *}*/
	Vector<double> tmp(2);
	switch(i%3){
		case 0:
			{
				/*{!therefore for these sites
				 *    f--g
				 *     \/ 
				 *     1  
				 *     /\
				 * 2--x--0
				 *   /
				 *  3 
				 *}*/
				switch(dir){
					case 0:
						{
							tmp(0) = this->dir_nn_LxLy_(0,0);
							tmp(1) = this->dir_nn_LxLy_(0,1);
						}break;
					case 1:
						{
							tmp(0) = this->dir_nn_LxLy_(2,0);
							tmp(1) = this->dir_nn_LxLy_(2,1);
						}break;
					case 2:
						{
							tmp(0) =-this->dir_nn_LxLy_(0,0);
							tmp(1) =-this->dir_nn_LxLy_(0,1);
						}break;
					case 3:
						{
							tmp(0) =-this->dir_nn_LxLy_(2,0);
							tmp(1) =-this->dir_nn_LxLy_(2,1);
						}break;
				}
			}break;
		case 1:
			{
				/*{!therefore for these sites
				 * f--g
				 *  \/ 
				 *  1  
				 *  /\
				 * 2--x--0
				 *     \
				 *      3 
				 *}*/
				switch(dir){
					case 0:
						{
							tmp(0) = this->dir_nn_LxLy_(0,0);
							tmp(1) = this->dir_nn_LxLy_(0,1);
						}break;
					case 1:
						{
							tmp(0) = this->dir_nn_LxLy_(1,0);
							tmp(1) = this->dir_nn_LxLy_(1,1);
						}break;
					case 2:
						{
							tmp(0) =-this->dir_nn_LxLy_(0,0);
							tmp(1) =-this->dir_nn_LxLy_(0,1);
						}break;
					case 3:
						{
							tmp(0) =-this->dir_nn_LxLy_(1,0);
							tmp(1) =-this->dir_nn_LxLy_(1,1);
						}break;
				}
			}break;
		case 2:
			{
				/*{!therefore for these sites
				 * 1--0
				 *  \/ 
				 *  x  
				 *  /\
				 * 2--3
				 *}*/
				switch(dir){
					case 0:
						{
							tmp(0) = this->dir_nn_LxLy_(2,0);
							tmp(1) = this->dir_nn_LxLy_(2,1);
						}break;
					case 1:
						{
							tmp(0) = this->dir_nn_LxLy_(1,0);
							tmp(1) = this->dir_nn_LxLy_(1,1);
						}break;
					case 2:
						{
							tmp(0) =-this->dir_nn_LxLy_(2,0);
							tmp(1) =-this->dir_nn_LxLy_(2,1);
						}break;
					case 3:
						{
							tmp(0) =-this->dir_nn_LxLy_(1,0);
							tmp(1) =-this->dir_nn_LxLy_(1,1);
						}break;
				}
			}break;
	}
	return tmp;
}

template<typename Type>
void Kagome<Type>::try_neighbourg(Vector<double>& tn, unsigned int const& i) const {
	if(i%this->xloop_){
		switch(i%3){
			case 0:
				{
					tn(0) += this->dir_nn_LxLy_(0,0)-this->dir_nn_LxLy_(1,0); 
					tn(1) += this->dir_nn_LxLy_(0,1)-this->dir_nn_LxLy_(1,1); 
				}break;
			case 1:
				{
					tn(0) += this->dir_nn_LxLy_(0,0);
					tn(1) += this->dir_nn_LxLy_(0,1);
				}break;
			case 2:
				{
					tn(0) += this->dir_nn_LxLy_(1,0);
					tn(1) += this->dir_nn_LxLy_(1,1);
				}break;
		}
	} else {
		tn(0) += this->dir_nn_LxLy_(2,0);
		tn(1) += this->dir_nn_LxLy_(2,1);
	}
}

template<typename Type>
Vector<double> Kagome<Type>::set_linear_jump() const {
	/*{!As the minimal unit cell contains 3 sites, the linear size is given by
	 * going from 0 to 3
	 *
	 *  2
	 *  /\
	 * 0--1--3
	 *}*/
	Vector<double> tmp(2);
	tmp(0) = 2.0;
	tmp(1) = 0;
	return tmp;
}
/*}*/
#endif
