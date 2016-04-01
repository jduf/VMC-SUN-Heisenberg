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
		void create_obs(unsigned int const& which_obs);
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Given N and m, save the best simulation in a text file for any n*/
		std::string extract_level_3();
};

/*{constructor*/
template<typename Type>
Ladder<Type>::Ladder(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,3,filename)
{
	/*!(*this) has been created via System(System const& s), this->J_ and
	 * this->links_ are already a copy of s. therefore if s has correctly
	 * defined J_ and links_, there is no need to recompute them for (*this).
	 * if s has undefined links_ and if J_ is of size 2, then this->links_ and
	 * this_->J_ should be computed*/
	if(this->status_==3){
		/*!create the links if necessary*/
		if(this->ref_(4)==2){
			Vector<unsigned int> tmp(2);
			tmp(0) = 2;
			tmp(1) = 1;
			this->create_energy_obs(tmp);
		} else {
			this->ref_(4) = 0;
			this->status_ = 2;
		}
		

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() == 2){
			Vector<double> tmp(this->J_);
			this->J_.set(this->obs_[0].nlinks());
			for (unsigned int i=0; i<this->J_.size();i++){
				if(i%3==1){ this->J_(i) = tmp(1); } //rungs (J⊥) -> sin(theta)
				else { this->J_(i) = tmp(0); }        //legs  (J‖) -> cos(theta)
			}
		}

		/*!fix the names for the bond energy*/
		if(this->J_.size()==this->obs_[0].nlinks()){
			std::string tmp("theta"+my::tostring(acos(this->J_(0))));
			this->filename_.replace(this->filename_.find("Juniform"),8,tmp);
			this->path_.replace(this->path_.find("Juniform"),8,tmp);
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : J_ has an incoherent size"<<std::endl; }
	} else {
		/*!if the ladder has a spuc_ equal to one, the creation is impossible
		 * but the construction should be silent to have a nice display when
		 * used with Analyse* */
		if(this->spuc_!=1){ std::cerr<<__PRETTY_FUNCTION__<<" creation is problematic"<<std::endl; }
	}
}

template<typename Type>
Ladder<Type>::~Ladder() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Ladder<Type>::create_obs(unsigned int const& which_obs){
	switch(which_obs){
		case 0:
			{ for(unsigned int i(1);i<3;i++){ create_obs(i); } }break;
		case 1:
			{
				unsigned int idx(this->obs_.size());
				this->obs_.push_back(Observable("Bond energy",1,this->z_*this->spuc_/2,this->obs_[0].nlinks()));
				this->obs_[idx].remove_links();
			}break;
		case 2:
			{
				unsigned int m(this->n_/2);
				unsigned int nval(this->n_/2);
				unsigned int nlinks(m*nval);
				unsigned int idx(this->obs_.size());
				this->obs_.push_back(Observable("S_10*S1i",2,nval,nlinks));
				this->obs_.push_back(Observable("S_10*S2i",2,nval,nlinks));
				this->obs_.push_back(Observable("S_20*S1i",2,nval,nlinks));
				this->obs_.push_back(Observable("S_10*S2i",2,nval,nlinks));
				for(unsigned int i(0);i<m;i++){
					for(unsigned int j(0);j<nval;j++){
						/*obs_[1]=S_10*S_1i*/
						this->obs_[idx](i*nval+j,0) = 2*i;
						this->obs_[idx](i*nval+j,1) = (2*(i+j))%this->n_;
						this->obs_[idx](i*nval+j,2) = j;
						/*obs_[2]=S_10*S_2i*/
						this->obs_[idx+1](i*nval+j,0) = 2*i;
						this->obs_[idx+1](i*nval+j,1) = (2*(i+j)+1)%this->n_;
						this->obs_[idx+1](i*nval+j,2) = j;
						/*obs_[3]=S_20*S_1i*/
						this->obs_[idx+2](i*nval+j,0) = 2*i+1;
						this->obs_[idx+2](i*nval+j,1) = (2*(i+j))%this->n_;
						this->obs_[idx+2](i*nval+j,2) = j;
						/*obs_[4]=S_20*S_2i*/
						this->obs_[idx+3](i*nval+j,0) = 2*i+1;
						this->obs_[idx+3](i*nval+j,1) = (2*(i+j)+1)%this->n_;
						this->obs_[idx+3](i*nval+j,2) = j;
					}
				}
			}break;
	}
}

template<typename Type>
Matrix<int> Ladder<Type>::get_neighbourg(unsigned int const& i) const {
	Matrix<int> nb(this->z_,3);
	if(i%2){//!odd number => upper part. 0:right, 1:down, 2:left
		if(i+1 != this->n_){
			nb(0,0) = i+2;
			nb(0,1) = 0;
			nb(0,2) = 0;
		} else {
			nb(0,0) = 1;
			nb(0,1) = 0;
			nb(0,2) = 1;
		}
		nb(1,0) = i-1;
		nb(1,1) = 1;
		nb(1,2) = 0;
		if(i != 1){
			nb(2,0) = i-2;
			nb(2,1) = 2;
			nb(2,2) = 0;
		} else {
			nb(2,0) = this->n_-1;
			nb(2,1) = 2;
			nb(2,2) = 1;
		}
	} else {//!even number => lower part. 0:right, 1:up, 2:left
		if(i+2 != this->n_){
			nb(0,0) = i+2;
			nb(0,1) = 0;
			nb(0,2) = 0;
		} else{
			nb(0,0) = 0;
			nb(0,1) = 0;
			nb(0,2) = 1;
		}
		nb(1,0) = i+1;
		nb(1,1) = 1;
		nb(1,2) = 0;
		if(i){
			nb(2,0) = i-2;
			nb(2,1) = 2;
			nb(2,2) = 0;
		} else {
			nb(2,0) = this->n_-2;
			nb(2,1) = 2;
			nb(2,2) = 1;
		}
	}
	return nb;
}

template<typename Type>
std::string Ladder<Type>::extract_level_3(){
	(*this->read_)>>this->obs_[0][0];
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->obs_[0][0]<<IOFiles::endl;

	return this->filename_;
}
/*}*/
#endif
