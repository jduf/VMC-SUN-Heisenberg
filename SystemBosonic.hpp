#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "MCSystem.hpp"
#include "Bosonic.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic : public Bosonic<Type>, public MCSystem<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemBosonic(Bosonic<Type> const& S, Rand& seed);
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
		/*!Forbids copy*/
		SystemBosonic(SystemBosonic const& S);
		/*!Forbids assignment*/
		SystemBosonic& operator=(SystemBosonic const& S);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic(Bosonic<Type> const& S, Rand& seed):
	System(S),
	Bosonic<Type>(S),
	MCSystem<Type>(S,seed)
{
	std::cerr<<"SystemBosonic will need to check everything"<<std::endl;
	Vector<unsigned int> available(this->n_);
	unsigned int N_as(this->n_);
	unsigned int site(0);
	for(unsigned int i(0); i < this->n_; i++){
		available(i) = i;
	}
	for(unsigned int c(0); c<this->N_; c++){
		for(unsigned int i(0); i < this->M_(c); i++){
			site = this->rnd_.get(N_as);
			this->s_(available(site),0) = c;
			this->s_(available(site),1) = c*this->M_(c)+i;
			for(unsigned int j(site); j+1 < N_as; j++){
				available(j) = available(j+1);
			}
			N_as--;
		}
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemBosonic<Type>::ratio(){
	if(this->new_c[0] == this->new_c[1]){
		/*!the minus sign is required because two particles are exchanged
		 *(Marshall-Peirels sign rule)*/
		/*not sure the -1 is correct*/
		return 1.0;
	} else {
		Type omegab_a(0.0);//omega_next/omega_current
		/*next state*/
		omegab_a = omega_(sl_(this->new_s[1]),this->new_c[0])
			* omega_(sl_(this->new_s[0]),this->new_c[1]); 
		/*current state*/
		omegab_a /= omega_(sl_(this->new_s[0]),this->new_c[0])
			* omega_(sl_(this->new_s[1]),this->new_c[1]);

		double jastrow(0.0);
		unsigned int c0,c1;
		for(unsigned int i(0);i<this->nn_.col();i++){
			c0=this->s_(this->nn_(this->new_s[0],i),0);
			c1=this->s_(this->nn_(this->new_s[1],i),0);
			if(nn_(this->new_s[0],i) != this->new_s[1]){
				jastrow += nu_(i, cc_(this->new_c[0], c0));
				jastrow -= nu_(i, cc_(this->new_c[1], c0));
			} else {
				jastrow += nu_(i, cc_(this->new_c[0], c0))/2.0;
				jastrow -= nu_(i, cc_(this->new_c[1], this->new_c[0]))/2.0;
			}
			if(nn_(this->new_s[1],i) != this->new_s[0]){
				jastrow += nu_(i, cc_(this->new_c[1], c1));
				jastrow -= nu_(i, cc_(this->new_c[0], c1));
			} else {
				jastrow += nu_(i, cc_(this->new_c[1], c1))/2.0;
				jastrow -= nu_(i, cc_(this->new_c[0], this->new_c[1]))/2.0;
			}
		}
		return exp(jastrow)*omegab_a;
	}
}
/*}*/
#endif
