#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "MCSystem.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic : public MCSystem<Type>, Bosonic<Type>{
	public:
		SystemBosonic(CreateSystem const& cs, unsigned int const& type);
		SystemBosonic();
		~SystemBosonic();

		//{Description
		/*!Computes the ratio of the two Jastrow factor related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio
		 */ //}
		Type ratio();

		void set(CreateSystem const& cs, unsigned int const& type);

	private:
		/*!Forbids copy constructor*/
		SystemBosonic(SystemBosonic const& S);
		/*!Forbids assignment operator*/
		SystemBosonic& operator=(SystemBosonic const& S);

		void init();
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic(CreateSystem const& cs, unsigned int const& type) {
	set(cs,type);
	std::cout<<"ok normal SystemBosonic"<<std::endl;
}
template<typename Type>
SystemBosonic<Type>::SystemBosonic(){
	std::cout<<"ok default SystemBosonic"<<std::endl;
}

template<typename Type>
void SystemBosonic<Type>::init(){
	std::cout<<"Bosonic"<<std::endl;
	Vector<unsigned int> available(this->n_);
	unsigned int N_as(this->n_);
	unsigned int site(0);
	for(unsigned int i(0); i < this->n_; i++){
		available(i) = i;
	}
	for(unsigned int c(0); c<this->N_; c++){
		for(unsigned int i(0); i < this->M_; i++){
			site = this->rnd_->get(N_as);
			this->s_(available(site),0) = c;
			this->s_(available(site),1) = c*this->M_+i;
			for(unsigned int j(site); j+1 < N_as; j++){
				available(j) = available(j+1);
			}
			N_as--;
		}
	}
}

template<typename Type>
SystemBosonic<Type>::~SystemBosonic(){}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void SystemBosonic<Type>::set(CreateSystem const& cs, unsigned int const& type){ 
	std::cout<<"SystemBosonic::set called"<<std::endl;
	/*init Bosonic*/
	this->nn_ = (cs.get<Type>())->get_nn();
	this->cc_ = (cs.get<Type>())->get_cc();
	this->nu_ = (cs.get<Type>())->get_nu();
	this->sl_ = (cs.get<Type>())->get_sl();
	this->omega_ = (cs.get<Type>())->get_omega();
	/*init MCSystem*/
	this->type_ = type;
	/*init System*/
	this->n_ = (cs.get<Type>())->get_n();
	this->N_ = (cs.get<Type>())->get_N();
	this->m_ = (cs.get<Type>())->get_m();
	this->M_ = (cs.get<Type>())->get_M();
	this->bc_ = (cs.get<Type>())->get_bc();
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
