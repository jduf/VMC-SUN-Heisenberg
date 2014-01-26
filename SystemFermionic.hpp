#ifndef DEF_SYSTEMFERMIONIC
#define DEF_SYSTEMFERMIONIC

#include "System.hpp"
#include "Lapack.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemFermionic : public System<Type>{
	public:
		/*!Creates a SystemFermionic without any parameters set*/
		SystemFermionic();
		/*!delete all the variables dynamically allocated*/
		~SystemFermionic();

		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - calls System<Type>::init()
		 * - sets tmp, U, Ut to the correct size 
		 * - creates an random initial state
		 * - if the sate is allowed, compute its related Ainv matrices
		*/ //}
		unsigned int init(Container const& input, unsigned int const& thread);

		/*!Call System<Type>::swap() and set row and new_ev*/
		void swap();

		/*!Calls System<Type>::swap(unsigned int const& s0, unsigned int const&
		 * s1) and set row and new_ev*/
		void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);

		//{Description
		/*!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 *
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio using
		 *   the determinant lemma
		 */ //}
		Type ratio();

		//{Description
		/*!Calls System<Type>::update() and then
		 *
		 * - updates the Ainv_ matrices
		 * - updates s_(new_s[i],1)
		 */ //}
		void update();

		void print();

	private:
		/*!Forbids copy constructor*/
		SystemFermionic(SystemFermionic const& S);
		/*!Forbids assignment operator*/
		SystemFermionic& operator=(SystemFermionic const& S);

		Matrix<unsigned int> row_;
		Matrix<Type> EVec_;		//!< det(A) <=> <GS|a>
		Matrix<Type> *Ainv_;	//!< inverse of A
		Matrix<Type> tmp;		//!< temporary matrix used during the update 
		Type w[2];				//!< det(W)= d = determinant ratios of <GS|a>/<GS|b> ; W=(w11,0,0,w22)
		Type d;					//!< d=Det(W) : determinant ratios of <GS|a>/<GS|b> 
		unsigned int new_r[2];	//!< rows of the Ainv_ matrix that are modified (the rows of the related A matrix are modified)
		unsigned int new_ev[2]; //!< new selected rows of the EVec matrix
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SystemFermionic<Type>::SystemFermionic():
	System<Type>(),
	Ainv_(NULL),
	d(0.0)
{}

template<typename Type>
SystemFermionic<Type>::~SystemFermionic(){
	delete[] Ainv_;
}

template<typename Type>
unsigned int SystemFermionic<Type>::init(Container const& input, unsigned int const& thread){
	System<Type>::init(input,thread);
	EVec_= input.get<Matrix<Type> >("EVec");
	Ainv_ = new Matrix<Type>[this->N_];
	for(unsigned int i(0); i < this->N_; i++){
		Ainv_[i].set(this->pps_*this->m_,this->pps_*this->m_);
	}
	tmp.set(this->pps_*this->m_,this->pps_*this->m_);
	row_.set(this->n_,this->pps_);

	Vector<int> ipiv;
	unsigned int TRY_MAX(100);
	unsigned int l(0);
	double rcn(0.0);
	do {
		unsigned int k(0);
		for(unsigned int i(0); i<this->n_; i++){
			for(unsigned int j(0); j<this->pps_; j++){
				this->s_(i,j) = ++k % this->N_;
			}
		}
		std::cout<<this->s_<<std::endl;
		for(unsigned int i(0); i<10; i++){
			swap();
			this->s_(this->new_s[0],this->new_p[0]) = this->new_c[1];
			this->s_(this->new_s[1],this->new_p[1]) = this->new_c[0];
		}
		std::cout<<"------------------"<<std::endl;
		std::cout<<this->s_<<std::endl;

		unsigned int c(0);
		unsigned int r;
		Vector<unsigned int> a(this->N_,0);
		for(unsigned int s(0); s < this->n_; s++){
			for(unsigned int p(0); p < this->pps_; p++){
				c = this->s_(s,p);
				r = c*this->n_+s;
				for(unsigned int j(0); j < this->pps_*this->m_; j++){
					Ainv_[c](a(c),j) = EVec_(r,j);
				}
				row_(s,p) = a(c);
				a(c)++;
			}
		}

		std::cout<<row_<<std::endl;
		std::cout<<Ainv_[2]<<std::endl;
		std::cout<<std::endl;
		std::cout<<EVec_<<std::endl;
		for(unsigned int c(0); c < this->N_; c++){
			Lapack<Type> inv(&Ainv_[c],false,'G');
			ipiv = inv.is_singular(rcn);
			if(!ipiv.ptr()){
				c = this->N_;
			} else {
				inv.inv(ipiv);
			}
		}

	} while (!ipiv.ptr() && ++l<TRY_MAX);

	if(l==TRY_MAX){
		std::cerr<<"sorry, the thread will not be lunched because no initial state was found"<<std::endl;
		return 0;
	} else {
		std::cerr<<"yeah ! initial state found"<<std::endl;
		return 1;
	}
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void SystemFermionic<Type>::update(){
	///*update the sites*/
	System<Type>::update();
	row_(this->new_s[0],this->new_p[0]) = new_r[1];
	row_(this->new_s[1],this->new_p[1]) = new_r[0];

	Type t_tmp;
	for(unsigned int c(0);c<2;c++){
		for(unsigned int j(0);j<this->pps_*this->m_;j++){
			if(new_r[c] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<this->pps_*this->m_;k++){
				t_tmp += EVec_(new_ev[c],k)*Ainv_[this->new_c[c]](k,j);
			}
			for(unsigned int i(0);i<this->pps_*this->m_;i++){
				tmp(i,j) = t_tmp*Ainv_[this->new_c[c]](i,new_r[c])/w[c];
			}
		}
		Ainv_[this->new_c[c]] -= tmp;
		//ev(this->new_c[c]*this->m_+new_r[c]) = new_ev[c];
	}
}

template<typename Type>
void SystemFermionic<Type>::swap(){
	System<Type>::swap();
	new_r[0] = row_(this->new_s[0],this->new_p[0]); 
	new_r[1] = row_(this->new_s[1],this->new_p[1]); 
	new_ev[0] = this->new_c[0]*this->n_ + this->new_s[1];
	new_ev[1] = this->new_c[1]*this->n_ + this->new_s[0];
}

template<typename Type>
void SystemFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	System<Type>::swap(s0,s1,p0,p1);
	new_r[0] = row_(s0,p0);
	new_r[1] = row_(s1,p1); 
	new_ev[0] = this->new_c[0]*this->n_ + s1;
	new_ev[1] = this->new_c[1]*this->n_ + s0;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemFermionic<Type>::ratio(){
	//std::cout<<new_ev[0]<<" "<<new_ev[1]<<std::endl;
	//std::cout<<new_r[0]<<" "<<new_r[1]<<std::endl;
	if(this->new_c[0] == this->new_c[1]){
		/*if swap() allow exchanges of the same color, then w and d has to be
		 * computed here*/
		return 1.0;
	} else {
		for(unsigned int i(0);i<2;i++){
			w[i] = 0.0;
			for(unsigned int k(0);k<this->pps_*this->m_;k++){
				w[i] += EVec_(new_ev[i],k)*Ainv_[this->new_c[i]](k,new_r[i]);
			}
		}
		d=w[0]*w[1];
		/*will simply need to multiply d by the jastrow*/
		if( std::abs(d) < 1e-10 ){ return 0.0; }
		else { return -d; }
	}
}

template<typename Type>
void SystemFermionic<Type>::print(){
	std::cout<<"=========================="<<std::endl;

	for(unsigned int i(0); i < this->n_; i++){
		std::cout<<i<<"("<<this->s_(i,0)<<","<<row_(i,0)<<") ";
	}
	std::cout<<std::endl;
	std::cout<<this->new_s[0]<<" "<<this->new_s[1]<<std::endl;
	std::cout<<"=========================="<<std::endl;
}
/*}*/
#endif
