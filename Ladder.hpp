#ifndef DEF_LADDER
#define DEF_LADDER

#include "System1D.hpp"

/*{Description*/
/*! 
 * This is our ladder.
 * 
 * 		1___3___5___7__...__n-1
 * 		|   |   |   |        |       
 * 		0___2___4___6__...__n-2
 */
/*}*/
template<typename Type>
class Ladder: public System1D<Type>{
	public:
		/*{Description*/
		/*!Constructor of a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(3,filename), to construct a system with 3 links
		 * per sites */
		/*}*/
		Ladder(unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Ladder()=0;

		Vector<double> compute_J(Vector<double> const& J) const;

	protected:
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Given N and m, save the best simulation in a text file for any n*/
		std::string extract_level_3();
};

template<typename Type>
Ladder<Type>::Ladder(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,3,filename)
{
	if(this->status_==2){ 
		Vector<unsigned int> l(2);
		l(0) = 2;
		l(1) = 1;
		this->compute_links(l);
	}
}

template<typename Type>
Ladder<Type>::~Ladder() = default;

template<typename Type>
Vector<double> Ladder<Type>::compute_J(Vector<double> const& J) const {
	Vector<double> tmp(this->links_.row());
	if(J.size() == 2){
		for (unsigned int i=0; i<this->links_.row() ; i++){
			if (i%3==1){ tmp(i) = J(0); } //rungs (J⊥)
			else{ tmp(i) = J(1); } //(J‖)
		}
	} else {
		std::cerr<<"Vector<double> const& create_J(Vector<double> const& J) : need J.size() == 2"<<std::endl;
	}
	return tmp;
}

template<typename Type>
Matrix<int> Ladder<Type>::get_neighbourg(unsigned int const& i) const {
	Matrix<int> nb(this->z_,2,1);
	if(i%2){// odd number => upper part. 0:right, 1:down, 2:left
		if(i+1 != this->n_){nb(0,0) = i+2;}
		else {
			nb(0,0) = 1;
			nb(0,1) = this ->bc_;
		}
		nb(1,0)=i-1;
		if(i != 1){nb(2,0)=i-2;}
		else {
			nb(2,0) = this-> n_-1;
			nb(2,1) = this-> bc_;
		}
	} else {// even number => lower part. 0:right, 1:up, 2:left
		if(i+2 != this->n_){nb(0,0) = i+2;}
		else {
			nb(0,0) = 0;
			nb(0,1) = this ->bc_;
		}
		nb(1,0)=i+1;
		if(i){nb(2,0)=i-2;}
		else {
			nb(2,0) = this-> n_-2;
			nb(2,1) = this-> bc_;
		}
	}
	return nb;
}

template<typename Type>
std::string Ladder<Type>::extract_level_3(){
	(*this->read_)>>this->E_;
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->E_<<IOFiles::endl;

	return this->filename_;
}
#endif
