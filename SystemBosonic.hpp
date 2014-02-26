#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic : public System<Type>{
	public:
		/*!Creates a SystemBosonic without any parameters set*/
		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - calls System<Type>::init()
		 * - creates an random initial state
		 */ //}
		SystemBosonic(CreateSystem const& CS, unsigned int const& thread);
		/*!delete all the variables dynamically allocated*/
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

		void print();

	private:
		/*!Forbids copy constructor*/
		SystemBosonic(SystemBosonic const& S);
		/*!Forbids assignment operator*/
		SystemBosonic& operator=(SystemBosonic const& S);

		Matrix<unsigned int> nn_; //!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;
		Matrix<double> nu_;
		Vector<unsigned int> sl_;
		Matrix<Type> omega_;
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic(CreateSystem const& CS, unsigned int const& thread):
	System<Type>(CS,thread)
	//nu_(CS.get_nu()),
	//nn_(CS.get_nn()),
	//cc_(CS.get_cc()),
	//sl_(CS.get_sl()),
	//omega_(CS.get_omega())
{
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
		for(unsigned int i(0);i<nn_.col();i++){
			c0=this->s_(nn_(this->new_s[0],i),0);
			c1=this->s_(nn_(this->new_s[1],i),0);
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

template<typename Type>
void SystemBosonic<Type>::print(){
	for(unsigned int i(0); i < sqrt(this->n_); i++){
		for(unsigned int j(0); j < sqrt(this->n_); j++){
			std::cout<<this->s_(j+i*int(sqrt(this->n_)),0)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"=========================="<<std::endl;
}
/*}*/
#endif
