#ifndef DEF_SYSTEMFERMIONIC
#define DEF_SYSTEMFERMIONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemFermionic : public System<Type>{
	public:
		/*!create a SystemFermionic without any parameters set*/
		SystemFermionic();
		/*!delete all the variables dynamically allocated*/
		~SystemFermionic();

		/*!exchanges two particles of different color */
		void swap();
		/*!exchanges particle on site s1 with the one on site s2*/
		void swap(unsigned int const& s0, unsigned int const& s1);
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
		/*!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio using
		 *   the determinant lemma
		 */
		//}
		Type ratio();
		//{Description
		/*!Updates the state if the condition given by the System::ratio()
		 * method is accepted. The update consists of :
		 *
		 * - computes the Ainv_ matrices
		 * - updates the configuration : s
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
		SystemFermionic(SystemFermionic const& S);
		/*!Forbids assignment operator*/
		SystemFermionic& operator=(SystemFermionic const& S);

		Rand* rnd;				//!< generator of random numbers 
		Matrix<Type> EVec_;		//!< det(A) <=> <GS|a>
		Matrix<Type> Ainv_;		//!< inverse of A
		Matrix<Type> U,Ut;		//!< temporary matrices used during the update 
		Matrix<Type> tmp;		//!< temporary matrix used during the update 
		Type w[4];				//!< W=(w11,w12,w21,w22)
		Type d;					//!< d=Det(W) : determinant ratios of <GS|a>/<GS|b> 
		unsigned int row[2];	//!< rows of the Ainv_ matrix that are modified (the rows of the related A matrix are modified)
		unsigned int new_ev[2]; //!< new selected rows of the EVec matrix

		/*could be deleted*/
		Vector<unsigned int> ev;//!< ev(color*m_+site)=row of EVec to select 
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemFermionic<Type>::SystemFermionic():
	System<Type>(),
	d(0.0)
{ }

template<typename Type>
SystemFermionic<Type>::~SystemFermionic(){ }
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
unsigned int SystemFermionic<Type>::init(Container const& input, unsigned int thread){
	System<Type>::init(input,thread);
	EVec_= input.get<Matrix<Type> >("EVec");
	Ainv_.set(this->n_,this->n_);
	tmp.set(this->n_,this->n_);
	Lapack<Type> inv(&Ainv_,false,'G');
	U.set(this->n_,2);
	Ut.set(2,this->n_);

	Vector<unsigned int> available(this->n_);
	Vector<int> ipiv;
	unsigned int N_as(this->n_);
	unsigned int site(0);
	unsigned int TRY_MAX(100);
	unsigned int l(0);
	double rcn(0.0);
	do {
		N_as = this->n_;
		for(unsigned int i(0); i < this->n_; i++){
			available(i) = i;
		}
		for(unsigned int c(0); c < this->N_; c++){
			for(unsigned int i(0); i < this->m_; i++){
				site = rnd->get(N_as);
				ev(c*this->m_+i) = c*this->n_+available(site);
				this->s_(available(site),0) = c;
				this->s_(available(site),1) = c*this->m_+i;
				for(unsigned int j(site); j+1 < N_as; j++){
					available(j) = available(j+1);
				}
				N_as--;
			}
		}
		for(unsigned int i(0); i < this->n_; i++){
			for(unsigned int j(0); j < this->n_; j++){
				Ainv_(i,j) = EVec_(ev(i),j);
			}
		}
		ipiv = inv.is_singular(rcn);
	} while (!ipiv.ptr() && ++l<TRY_MAX);
	
	if(l==TRY_MAX){
		std::cerr<<"sorry, the thread will not be lunched because no initial state was found"<<std::endl;
		return 0;
	} else {
		std::cerr<<"yeah ! initial state found"<<std::endl;
		inv.inv(ipiv);
		return 1;
	}
}

template<typename Type>
void SystemFermionic<Type>::update(){
	///*update the sites*/
	System<Type>::update();
	s_(this->new_s[0],1) = row[1];
	s_(this->new_s[1],1) = row[0];

	/*compute the inverse of W*/
	Type t;
	t = w[0];
	w[0] = w[3]/d;
	w[3] = t/d;
	w[1] = -w[1]/d;
	w[2] = -w[2]/d;

	/*apply the woodbury formula*/
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<2;j++){
			U(i,j) = 0.0;
			for(unsigned int k(0);k<2;k++){
				U(i,j) += Ainv_(i,row[k])*w[k+j*2]; // Ainv.U.Winv
			}
		}
	}
	for(unsigned int i(0);i<2;i++){
		for(unsigned int j(0);j<this->n_;j++){
			Ut(i,j) = 0.0;
			for(unsigned int k(0);k<this->n_;k++){
				Ut(i,j) += EVec_(new_ev[i],k)*Ainv_(k,j); // Ut.Ã.Ainv
			}
		}
	}
	Ut(0,row[0]) -= 1.0;
	Ut(1,row[1]) -= 1.0;
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			tmp(i,j) = 0.0;
			for(unsigned int k(0);k<2;k++){
				tmp(i,j) += U(i,k)*Ut(k,j);
			}
		}
	}
	Ainv_ -= tmp;

	/*is only useful for SystemFermionic::print() : */
	ev(row[0]) = new_ev[0];
	ev(row[1]) = new_ev[1];
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
void SystemFermionic<Type>::swap(){
	System<Type>::swap();
	row[0] = this->s_(this->new_s[0],1); 
	row[1] = this->s_(this->new_s[1],1); 
	new_ev[0] = this->color[0]*this->n_ + this->new_s[1];
	new_ev[1] = this->color[1]*this->n_ + this->new_s[0];
}

template<typename Type>
void SystemFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1){
	System<Type>::swap(s0,s1);
	row[0] = this->s_(s0,1);
	row[1] = this->s_(s1,1); 
	new_ev[0] = this->color[0]*this->n_ + s1;
	new_ev[1] = this->color[1]*this->n_ + s0;
}

template<typename Type> 
Type SystemFermionic<Type>::ratio(){
	/*if swap() allow exchanges of the same color, then w and d has to be
	 * computed here*/
	if(this->color[0] == this->color[1]){
		/*!the minus sign is required because two particles are exchanged
		 *(Marshall-Peirels sign rule)*/
		return -1.0;
	} else {
		for(unsigned int i(0);i<2;i++){
			for(unsigned int j(0);j<2;j++){
				w[i+2*j] = 0.0;
				for(unsigned int k(0);k<this->n_;k++){
					w[i+j*2] += EVec_(new_ev[i],k) * Ainv_(k,row[j]);
				}
			}
		}
		d = w[0]*w[3]-w[1]*w[2];
		if( std::abs(d) < 1e-14 ){ return 0.0; }
		else { return d; }
	}
}

template<typename Type>
void SystemFermionic<Type>::print(){
	std::cout<<"=========================="<<std::endl;
	Matrix<Type> A(this->n_,this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			A(i,j) = EVec_(ev(i),j);
		}
	}
	std::cout<<"inv "<<(A*Ainv_).diag().chop()<<std::endl;
	//std::cout<<"ev :"<<ev<<std::endl;
	//std::cout<<"s : ";
	//for(unsigned int i(0); i < this->n_; i++){
	//std::cout<<i<<"("<<s(i,0)<<","<<s(i,1)<<") ";
	//}
	std::cout<<"=========================="<<std::endl;
}
/*}*/
#endif
