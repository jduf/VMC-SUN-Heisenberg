#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "System1D.hpp"

template<typename Type>
class Chain: public System1D<Type>{
	public:
		/*{Description*/
		/*!Constructor of a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(2,filename), to construct a system with 2 links
		 * per sites */
		/*}*/
		Chain(unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Chain()=0;

	protected:
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int i) const;
		std::string extract_level_3();
};

template<typename Type>
Chain<Type>::Chain(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,2,filename)
{
	if(this->status_==1){ this->compute_links(); }
}

template<typename Type>
Chain<Type>::~Chain(){}

template<typename Type>
Matrix<int> Chain<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	if( i != this->n_-1){ nb(0,0) = i+1;}
	else {
		nb(0,0) = 0;
		nb(0,1) = this->bc_;
	}
	if( i != 0){ nb(1,0) = i-1;}
	else {
		nb(1,0) = this->n_-1;
		nb(1,1) = this->bc_;
	}
	return nb;
}

template<typename Type>
std::string Chain<Type>::extract_level_3(){
	double polymerization_strength;
	(*this->read_)>>this->E_>>polymerization_strength;
	(*this->data_write_)<<this->n_<<" "<<this->E_<<" "<<polymerization_strength<<" "<<this->bc_<<IOFiles::endl;

	return this->filename_;
}
#endif
