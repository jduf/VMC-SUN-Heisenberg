#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "MCSystem.hpp"
#include "Bosonic.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic : public Bosonic<Type>, public MCSystem<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemBosonic(Bosonic<Type> const& S);
		/*!Destructor*/
		~SystemBosonic(){};

		//{Description
		/*!Computes the ratio of the two Jastrow factor related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio
		 */ //}
		Type ratio();

	private:
		/*!Forbids default*/
		SystemBosonic();
		/*!Forbids copy*/
		SystemBosonic(SystemBosonic const& S);
		/*!Forbids assignment*/
		SystemBosonic& operator=(SystemBosonic const& S);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic(Bosonic<Type> const& S):
	System(S),
	Bosonic<Type>(S),
	MCSystem<Type>(S)
{
	std::cerr<<"SystemBosonic<Type>::SystemBosonic(Bosonic<Type> const& S) : check everything"<<std::endl;
	std::cerr<<"works only for m=1"<<std::endl;
	this->status_--;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemBosonic<Type>::ratio(){
	if(this->new_c_[0] == this->new_c_[1]){
		/*!the minus sign is required because two particles are exchanged
		 *(Marshall-Peirels sign rule)*/
		/*not sure the -1 is correct*/
		return 1.0;
	} else {
		Type omegab_a(0.0);//this->omega_next/this->omega_current
		/*next state*/
		omegab_a = this->omega_(this->sl_(this->new_s_[1]),this->new_c_[0])
			* this->omega_(this->sl_(this->new_s_[0]),this->new_c_[1]); 
		/*current state*/
		omegab_a /= this->omega_(this->sl_(this->new_s_[0]),this->new_c_[0])
			* this->omega_(this->sl_(this->new_s_[1]),this->new_c_[1]);

		double jastrow(0.0);
		unsigned int c0,c1;
		for(unsigned int i(0);i<this->nn_.col();i++){
			c0=this->s_(this->nn_(this->new_s_[0],i),0);
			c1=this->s_(this->nn_(this->new_s_[1],i),0);
			if(this->nn_(this->new_s_[0],i) != this->new_s_[1]){
				jastrow += this->nu_(i, this->cc_(this->new_c_[0], c0));
				jastrow -= this->nu_(i, this->cc_(this->new_c_[1], c0));
			} else {
				jastrow += this->nu_(i, this->cc_(this->new_c_[0], c0))/2.0;
				jastrow -= this->nu_(i, this->cc_(this->new_c_[1], this->new_c_[0]))/2.0;
			}
			if(this->nn_(this->new_s_[1],i) != this->new_s_[0]){
				jastrow += this->nu_(i, this->cc_(this->new_c_[1], c1));
				jastrow -= this->nu_(i, this->cc_(this->new_c_[0], c1));
			} else {
				jastrow += this->nu_(i, this->cc_(this->new_c_[1], c1))/2.0;
				jastrow -= this->nu_(i, this->cc_(this->new_c_[0], this->new_c_[1]))/2.0;
			}
		}
		return exp(jastrow)*omegab_a;
	}
}
/*}*/
#endif
