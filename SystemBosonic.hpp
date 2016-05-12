#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "MCSystem.hpp"
#include "Bosonic.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic: public MCSystem, public Bosonic<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemBosonic(Bosonic<Type> const& S);
		/*!Constructor that reads from file*/
		SystemBosonic(IOFiles& r);
		/*!Default destructor*/
		~SystemBosonic() = default;
		/*{Forbidden*/
		SystemBosonic() = delete;
		SystemBosonic(SystemBosonic<Type>&) = delete;
		SystemBosonic& operator=(SystemBosonic<Type>) = delete;
		/*}*/

		/*{Description
		 * !Computes the ratio of the two Jastrow factor related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged a minus sign arises
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio
		 }*/
		double ratio();

		/*!Returns a copy of this instance*/
		std::unique_ptr<MCSystem> clone() const;
		/*!Set to most of the matrices to NULL*/
		void free_memory();

		void write(IOFiles& w) const;

	private:
		/*!Autorize copy only via clone()*/
		SystemBosonic(SystemBosonic<Type> const& S);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic(Bosonic<Type> const& S):
	System(S),
	MCSystem(S),
	Bosonic<Type>(S)
{
	std::cerr<<__PRETTY_FUNCTION__<<" : check everything, may only work for m=1"<<std::endl;
	status_=0;
}

template<typename Type>
SystemBosonic<Type>:: SystemBosonic(IOFiles& r):
	System(r),
	MCSystem(r),
	Bosonic<Type>(r)
{}

template<typename Type>
SystemBosonic<Type>::SystemBosonic(SystemBosonic<Type> const& S):
	System(S),
	MCSystem(S),
	Bosonic<Type>(S)
{}

template<typename Type>
std::unique_ptr<MCSystem> SystemBosonic<Type>::clone() const {
	return std::unique_ptr<SystemBosonic<Type> >(new SystemBosonic<Type>(*this));
}
/*}*/

/*void methods*/
/*{*/
template<typename Type>
void SystemBosonic<Type>::write(IOFiles& w) const {
	System::write(w);
	MCSystem::write(w);
	w<<this->sl_<<this->nn_<<this->cc_<<this->nu_<<this->omega_;
}

template<typename Type>
void SystemBosonic<Type>::free_memory(){
	std::cerr<<__PRETTY_FUNCTION__<<" : nothing to do"<<std::endl;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
double SystemBosonic<Type>::ratio(){
	if(new_c_[0] == new_c_[1]){
		/*!the minus sign is required because two particles are exchanged
		 *(Marshall-Peirels sign rule)*/
		/*not sure the -1 is correct*/
		return 1.0;
	} else {
		Type omegab_a(0.0);//this->omega_next/this->omega_current
		/*next state*/
		omegab_a = this->omega_(this->sl_(new_s_[1]),new_c_[0])
			* this->omega_(this->sl_(new_s_[0]),new_c_[1]);
		/*current state*/
		omegab_a /= this->omega_(this->sl_(new_s_[0]),new_c_[0])
			* this->omega_(this->sl_(new_s_[1]),new_c_[1]);

		double jastrow(0.0);
		unsigned int c0,c1;
		for(unsigned int i(0);i<this->nn_.col();i++){
			c0 = s_(this->nn_(new_s_[0],i),0);
			c1 = s_(this->nn_(new_s_[1],i),0);
			if(this->nn_(new_s_[0],i) != new_s_[1]){
				jastrow += this->nu_(i, this->cc_(new_c_[0], c0));
				jastrow -= this->nu_(i, this->cc_(new_c_[1], c0));
			} else {
				jastrow += this->nu_(i, this->cc_(new_c_[0], c0))/2.0;
				jastrow -= this->nu_(i, this->cc_(new_c_[1], new_c_[0]))/2.0;
			}
			if(this->nn_(new_s_[1],i) != new_s_[0]){
				jastrow += this->nu_(i, this->cc_(new_c_[1], c1));
				jastrow -= this->nu_(i, this->cc_(new_c_[0], c1));
			} else {
				jastrow += this->nu_(i, this->cc_(new_c_[1], c1))/2.0;
				jastrow -= this->nu_(i, this->cc_(new_c_[0], new_c_[1]))/2.0;
			}
		}
		return my::real(exp(jastrow)*omegab_a);
	}
}
/*}*/
#endif
