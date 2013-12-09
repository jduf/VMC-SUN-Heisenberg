#ifndef DEF_SYSTEMBOSONIC
#define DEF_SYSTEMBOSONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemBosonic : public System<Type>{
	public:
		/*!create a System without any parameters set*/
		SystemBosonic();
		/*!delete all the variables dynamically allocated*/
		~SystemBosonic();

		//{Description
		/*! This method creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - sets N, m, n, sts_ 
		 * - initialize the random number generator
		 * - creates an random initial state and computes
		 */
		//}
		unsigned int init(Container const& input, unsigned int thread);
		//{Description
		/*!Computes the ratio of the two Jastrow factor related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio
		 */
		//}
		Type ratio();
		//{Description
		/*!Updates the state if the condition given by the System::ratio()
		 * method is accepted. The update consists of :
		 *
		 * - updates the configuration : s
		 * - update the number of site that have different color
		 */
		//}
		void update();
		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		void measure(double& E_config);
		void print();

	private:
		/*!Forbids copy constructor*/
		SystemBosonic(SystemBosonic const& S);
		/*!Forbids assignment operator*/
		SystemBosonic& operator=(SystemBosonic const& S);

		Vector<unsigned int> fb_;//!< number of neighbour with a different color
		Matrix<unsigned int> are_neighbours_;
		double nu_;
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemBosonic<Type>::SystemBosonic():System<Type>(){std::cout<<"const ok"<<std::endl;}

template<typename Type>
SystemBosonic<Type>::~SystemBosonic(){
	print();
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int SystemBosonic<Type>::init(Container const& input, unsigned int thread){
	System<Type>::init(input,thread);
	fb_.set(this->n_,0);
	input.get("nu",nu_);
	are_neighbours_.set(this->n_,this->n_,0);

	Vector<unsigned int> available(this->n_);
	unsigned int N_as(this->n_);
	unsigned int site(0);
	for(unsigned int i(0); i < this->n_; i++){
		available(i) = i;
	}
	for(unsigned int color(0); color < this->N_; color++){
		for(unsigned int i(0); i < this->m_; i++){
			site = this->rnd->get(N_as);
			this->s_(available(site),0) = color;
			this->s_(available(site),1) = color*this->m_+i;
			for(unsigned int j(site); j+1 < N_as; j++){
				available(j) = available(j+1);
			}
			N_as--;
		}
	}

	for(unsigned int i(0); i < this->sts_.row(); i++){
		if(this->s_(this->sts_(i,0),0) != this->s_(this->sts_(i,1),0)){
			fb_(this->sts_(i,0))++;
			fb_(this->sts_(i,1))++;
		}
		are_neighbours_(this->sts_(i,0),this->sts_(i,1)) = 1;
		are_neighbours_(this->sts_(i,1),this->sts_(i,0)) = 1;
	}

	print();
	return 1;
}

template<typename Type>
void SystemBosonic<Type>::update(){
	///*update the sites*/
	System<Type>::update();

	/*update the jastrow*/
	if(this->color[0] != this->color[1]){
		for(unsigned int i(0); i<this->n_;i++){
			if(are_neighbours_(this->new_s[0],i)){
				if(this->s_(i,0) != this->color[1]) { fb_(i)++;}
				else { fb_(i)--;}
			}
			if(are_neighbours_(this->new_s[1],i)){
				if(this->s_(i,0) != this->color[0]) { fb_(i)++;}
				else { fb_(i)--;}
			}
		}
		fb_(this->new_s[0]) = 0;
		fb_(this->new_s[1]) = 0;
		for(unsigned int i(0); i<this->n_;i++){
			if(are_neighbours_(this->new_s[0],i) && this->s_(i,0) != this->color[1]) { fb_(this->new_s[0])++; }
			if(are_neighbours_(this->new_s[1],i) && this->s_(i,0) != this->color[0]) { fb_(this->new_s[1])++; }
		}
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemBosonic<Type>::ratio(){
	/*if swap() allow exchanges of the same color, then w and d has to be
	 * computed here*/
	if(this->color[0] == this->color[1]){
		/*!the minus sign is required because two particles are exchanged
		 *(Marshall-Peirels sign rule)*/
		return -1.0;
	} else {
		if(are_neighbours_(this->new_s[0],this->new_s[1])){
			return exp(nu_*(5.0-fb_(this->new_s[0])-fb_(this->new_s[1])));
		} else {
			return exp(nu_*(4.0-fb_(this->new_s[0])-fb_(this->new_s[1])));
		}
	}
}

template<typename Type>
void SystemBosonic<Type>::print(){
	for(unsigned int i(0); i < sqrt(this->n_); i++){
		for(unsigned int j(0); j < sqrt(this->n_); j++){
			std::cout<<"("<<this->s_(j+i*int(sqrt(this->n_)),0)<<","<<fb_(j+i*int(sqrt(this->n_)))<<") ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"=========================="<<std::endl;
}
/*}*/
#endif
