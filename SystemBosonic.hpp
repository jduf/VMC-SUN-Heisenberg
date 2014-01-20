#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic : public System<Type>{
	public:
		/*!Creates a SystemBosonic without any parameters set*/
		SystemBosonic();
		/*!delete all the variables dynamically allocated*/
		~SystemBosonic();

		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - calls System<Type>::init()
		 * - creates an random initial state
		 */ //}
		unsigned int init(Container const& input, unsigned int const& thread);

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
		void correlation(Matrix<unsigned int>* corr);

	private:
		/*!Forbids copy constructor*/
		SystemBosonic(SystemBosonic const& S);
		/*!Forbids assignment operator*/
		SystemBosonic& operator=(SystemBosonic const& S);

		Matrix<unsigned int> nn_; //!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;
		Matrix<double> nu_;
		Matrix<Type> omega_;
		Vector<unsigned int> sl_;
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic():
	System<Type>()
{}

template<typename Type>
SystemBosonic<Type>::~SystemBosonic(){}

template<typename Type>
unsigned int SystemBosonic<Type>::init(Container const& input, unsigned int const& thread){
	System<Type>::init(input,thread);
	nu_ = input.get<Matrix<double> >("nu");
	nn_ = input.get<Matrix<unsigned int> >("nn");
	cc_ = input.get<Matrix<unsigned int> >("cc");
	sl_ = input.get<Vector<unsigned int> >("sl");
	omega_ = input.get<Matrix<Type> >("omega");

	Vector<unsigned int> available(this->n_);
	unsigned int N_as(this->n_);
	unsigned int site(0);
	for(unsigned int i(0); i < this->n_; i++){
		available(i) = i;
	}
	for(unsigned int c(0); c<this->N_; c++){
		for(unsigned int i(0); i < this->m_; i++){
			site = this->rnd->get(N_as);
			this->s_(available(site),0) = c;
			this->s_(available(site),1) = c*this->m_+i;
			for(unsigned int j(site); j+1 < N_as; j++){
				available(j) = available(j+1);
			}
			N_as--;
		}
	}
	return 1;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemBosonic<Type>::ratio(){
	if(this->color[0] == this->color[1]){
		/*!the minus sign is required because two particles are exchanged
		 *(Marshall-Peirels sign rule)*/
		/*not sure the -1 is correct*/
		return 1.0;
	} else {
		Type omegab_a(0.0);//omega_next/omega_current
		/*next state*/
		omegab_a = omega_(sl_(this->new_s[1]),this->color[0])
			* omega_(sl_(this->new_s[0]),this->color[1]); 
		/*current state*/
		omegab_a /= omega_(sl_(this->new_s[0]),this->color[0])
			* omega_(sl_(this->new_s[1]),this->color[1]);

		double jastrow(0.0);
		unsigned int c0,c1;
		for(unsigned int i(0);i<nn_.col();i++){
			c0=this->s_(nn_(this->new_s[0],i),0);
			c1=this->s_(nn_(this->new_s[1],i),0);
			if(nn_(this->new_s[0],i) != this->new_s[1]){
				jastrow += nu_(i, cc_(this->color[0], c0));
				jastrow -= nu_(i, cc_(this->color[1], c0));
			} else {
				jastrow += nu_(i, cc_(this->color[0], c0))/2.0;
				jastrow -= nu_(i, cc_(this->color[1], this->color[0]))/2.0;
			}
			if(nn_(this->new_s[1],i) != this->new_s[0]){
				jastrow += nu_(i, cc_(this->color[1], c1));
				jastrow -= nu_(i, cc_(this->color[0], c1));
			} else {
				jastrow += nu_(i, cc_(this->color[1], c1))/2.0;
				jastrow -= nu_(i, cc_(this->color[0], this->color[1]))/2.0;
			}
		}
		return exp(jastrow)*omegab_a;
	}
}

template<typename Type> 
void SystemBosonic<Type>::correlation(Matrix<unsigned int>* corr){
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<nn_.col();j++){
			if(this->s_(i,0) == this->s_(nn_(i,j),0)){
				(*corr)(i,j)++;
			}
		}
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
