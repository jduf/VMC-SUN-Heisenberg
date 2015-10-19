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

	protected:
		void set_observables(int nobs);
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Given N and m, save the best simulation in a text file for any n*/
		std::string extract_level_3();
};

template<typename Type>
Ladder<Type>::Ladder(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,3,filename)
{
	/*!(*this) has been created via System(System const& s), this->J_ and
	 * this->links_ are already a copy of s. therefore if s has correctly
	 * defined J_ and links_, there is no need to recompute them for (*this).
	 * if s has undefined links_ and if J_ is of size 2, then this->links_ and
	 * this_->J_ should be computed*/
	if(this->status_==2){
		/*!create the links if necessary*/
		if(!this->obs_.size()){
			Vector<unsigned int> l(2);
			l(0) = 2;
			l(1) = 1;
			this->set_nn_links(l);
		}

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() == 2){
			Vector<double> tmp(this->J_);
			this->J_.set(this->obs_[0].nlinks());
			for (unsigned int i=0; i<this->J_.size();i++){
				if (i%3==1){ this->J_(i) = tmp(1); } //rungs (J⊥) -> sin(theta)
				else{ this->J_(i) = tmp(0); }        //legs  (J‖) -> cos(theta)
			}
		}

		/*!fix the names for the bond energy*/
		if(this->J_.size()==this->obs_[0].nlinks()){
			std::string tmp("theta"+my::tostring(acos(this->J_(0))));
			this->filename_.replace(this->filename_.find("Juniform"),8,tmp);
			this->path_.replace(this->path_.find("Juniform"),8,tmp);
		} else {
			std::cerr<<__PRETTY_FUNCTION__<<" : J_ has an incoherent size"<<std::endl;
		}
	} else {
		/*!if the ladder has a spuc_ equal to one, the creation is impossible
		 * but the construction should be silent to have a nice display when
		 * used with Analyse* */
		if(this->spuc_!=1){
			std::cerr<<__PRETTY_FUNCTION__<<" creation is problematic"<<std::endl;
		}
	}
}

template<typename Type>
Ladder<Type>::~Ladder() = default;

template<typename Type>
void Ladder<Type>::set_observables(int nobs){
	this->E_.set(50,5,false);

	if(nobs<0){ nobs = 6; }
	unsigned int nlinks;
	unsigned int nval;
	unsigned int m;
	if(nobs>0){/*bond energy*/
		nlinks = this->obs_[0].nlinks();
		nval = 3*this->spuc_/2;
		this->obs_[0].set(nval,50,5,false);
		for(unsigned int i(0);i<nlinks;i++){
			this->obs_[0](i,2) = i%nval;
		}
	}
	if(nobs>1){/*long range correlation*/
		m = this->spuc_/2;
		nval = this->n_;
		nlinks = m*nval;
		this->obs_.push_back(Observable(nlinks,nval,50,5,false));
		for(unsigned int i(0);i<m;i++){
			for(unsigned int j(0);j<nval;j++){
				this->obs_[1](i*nval+j,0) = (2*i)%this->n_;
				this->obs_[1](i*nval+j,1) = (2*i+j)%this->n_;
				this->obs_[1](i*nval+j,2) = j;
			}
		}
	}
	if(nobs==6){/*(anti)symmetric correlation*/
		m = this->n_/2;
		nval = this->n_/2;
		nlinks = m*nval;
		this->obs_.push_back(Observable(nlinks,nval,50,5,false));
		this->obs_.push_back(Observable(nlinks,nval,50,5,false));
		this->obs_.push_back(Observable(nlinks,nval,50,5,false));
		this->obs_.push_back(Observable(nlinks,nval,50,5,false));
		for(unsigned int i(0);i<m;i++){
			for(unsigned int j(0);j<nval;j++){
				//std::cout<<i<<" "<<j<<std::flush;
				/*obs_[2]=S_10*S_1i*/
				this->obs_[2](i*nval+j,0) = 2*i;
				this->obs_[2](i*nval+j,1) = (2*(i+j))%this->n_;
				this->obs_[2](i*nval+j,2) = j;
				/*obs_[3]=S_10*S_2i*/
				this->obs_[3](i*nval+j,0) = 2*i;
				this->obs_[3](i*nval+j,1) = (2*(i+j)+1)%this->n_;
				this->obs_[3](i*nval+j,2) = j;
				/*obs_[4]=S_20*S_1i*/
				this->obs_[4](i*nval+j,0) = 2*i+1;
				this->obs_[4](i*nval+j,1) = (2*(i+j))%this->n_;
				this->obs_[4](i*nval+j,2) = j;
				/*obs_[5]=S_20*S_2i*/
				this->obs_[5](i*nval+j,0) = 2*i+1;
				this->obs_[5](i*nval+j,1) = (2*(i+j)+1)%this->n_;
				this->obs_[5](i*nval+j,2) = j;
				//std::cout<<"| "<<this->obs_[2](i*m+j,0)<<" "<< this->obs_[2](i*m+j,1)<<" " <<this->obs_[2](i*m+j,2)<<std::endl;
				//std::cout<<"| "<<this->obs_[3](i*nlinks/2+j,0)<<" "<< this->obs_[3](i*nlinks/2+j,1)<<" " <<this->obs_[3](i*nlinks/2+j,2)<<std::endl;
				//std::cout<<"| "<<this->obs_[4](i*nlinks/2+j,0)<<" "<< this->obs_[4](i*nlinks/2+j,1)<<" " <<this->obs_[4](i*nlinks/2+j,2)<<std::endl;
				//std::cout<<"| "<<this->obs_[5](i*nlinks/2+j,0)<<" "<< this->obs_[5](i*nlinks/2+j,1)<<" " <<this->obs_[5](i*nlinks/2+j,2)<<std::endl;
			}
		}
	}
}

template<typename Type>
Matrix<int> Ladder<Type>::get_neighbourg(unsigned int const& i) const {
	Matrix<int> nb(this->z_,2,1);
	if(i%2){//!odd number => upper part. 0:right, 1:down, 2:left
		if(i+1 != this->n_){ nb(0,0) = i+2; }
		else{
			nb(0,0) = 1;
			nb(0,1) = this ->bc_;
		}
		nb(1,0)=i-1;
		if(i != 1){ nb(2,0)=i-2; }
		else{
			nb(2,0) = this-> n_-1;
			nb(2,1) = this-> bc_;
		}
	} else {//!even number => lower part. 0:right, 1:up, 2:left
		if(i+2 != this->n_){ nb(0,0) = i+2; }
		else{
			nb(0,0) = 0;
			nb(0,1) = this ->bc_;
		}
		nb(1,0)=i+1;
		if(i){ nb(2,0)=i-2; }
		else{
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
