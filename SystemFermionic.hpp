#ifndef DEF_SYSTEMFERMIONIC
#define DEF_SYSTEMFERMIONIC

#include "MCSystem.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class SystemFermionic : public MCSystem<Type>, Fermionic<Type>{
	public:
		//SystemFermionic(System* S,unsigned int const& thread, unsigned int const& type);
		SystemFermionic(System* S, Fermionic<Type>* bosonic);
		/*!Delete all the variables dynamically allocated*/
		~SystemFermionic();

		/*!Set row and new_ev*/
		void swap();
		/*!Set row and new_ev*/
		void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);

		//{Description
		/*!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio using
		 *   the determinant lemma
		 */ //}
		Type ratio();
		//{Description
		/*!Calls System::update() and then
		 * - updates the Ainv_ matrices
		 * - updates s_(new_s[i],1)
		 */ //}
		void update();

		void init();

	private:
		/*!Forbids copy constructor*/
		SystemFermionic(SystemFermionic const& S);
		/*!Forbids assignment operator*/
		SystemFermionic& operator=(SystemFermionic const& S);

		Matrix<unsigned int> row_;
		Matrix<Type> *Ainv_;	//!< inverse of A
		Matrix<Type> tmp_;		//!< temporary matrix used during the update 
		Type w[2];				//!< det(W)= d = determinant ratios of <GS|a>/<GS|b> ; W=(w11,0,0,w22)
		unsigned int new_r[2];	//!< rows of the Ainv_ matrix that are modified (the rows of the related A matrix are modified)
		unsigned int new_ev[2]; //!< new selected rows of the EVec matrix
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
SystemFermionic<Type>::SystemFermionic(System* S, Fermionic<Type>* fermionic):
	MCSystem<Type>(S),
	Fermionic<Type>(*fermionic),
	Ainv_(NULL)
{
	std::cout<<"ok SystemFermionic"<<std::endl;
}

template<typename Type>
void SystemFermionic<Type>::init(){
	Ainv_ = new Matrix<Type>[this->N_];
	for(unsigned int i(0); i < this->N_; i++){
		Ainv_[i].set(this->M_,this->M_);
	}
	tmp_.set(this->M_,this->M_);
	row_.set(this->n_,this->m_);

	Vector<int> ipiv;
	unsigned int TRY_MAX(100);
	unsigned int l(0);
	double rcn(0.0);
	unsigned int k(0);
	for(unsigned int i(0); i<this->n_; i++){
		for(unsigned int j(0); j<this->m_; j++){
			this->s_(i,j) = ++k % this->N_;
		}
	}
	do {
		for(unsigned int i(0); i<1000; i++){
			swap();
			this->s_(this->new_s[0],this->new_p[0]) = this->new_c[1];
			this->s_(this->new_s[1],this->new_p[1]) = this->new_c[0];
		}
		unsigned int c(0);
		unsigned int r;
		Vector<unsigned int> a(this->N_,0);
		for(unsigned int s(0); s < this->n_; s++){
			for(unsigned int p(0); p < this->m_; p++){
				c = this->s_(s,p);
				r = c*this->n_+s;
				for(unsigned int j(0); j < this->M_; j++){
					Ainv_[c](a(c),j) = this->EVec_(r,j);
				}
				row_(s,p) = a(c);
				a(c)++;
			}
		}
		for(unsigned int c(0); c < this->N_; c++){
			Lapack<Type> inv(&Ainv_[c],false,'G');
			ipiv = inv.is_singular(rcn);
			if(ipiv.ptr()){ inv.inv(ipiv); }
			else { c = this->N_; }
		}
	} while (!ipiv.ptr() && ++l<TRY_MAX);
	if(l!=TRY_MAX){
		this->ready_=true; /*2nd step successful*/
	} else {
		std::cerr<<"No initial state found after "<<TRY_MAX<<"trials"<<std::endl;
		this->ready_=false;
	}
	std::cout<<"system fermionic init ok"<<std::endl;
	std::cerr<<"need to initialize only if non degenerate"<<std::endl;
}

template<typename Type>
SystemFermionic<Type>::~SystemFermionic(){
	delete[] Ainv_;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void SystemFermionic<Type>::update(){
	MCSystem<Type>::update();
	row_(this->new_s[0],this->new_p[0]) = new_r[1];
	row_(this->new_s[1],this->new_p[1]) = new_r[0];

	Type t_tmp;
	for(unsigned int c(0);c<2;c++){
		for(unsigned int j(0);j<this->M_;j++){
			if(new_r[c] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<this->M_;k++){
				t_tmp += this->EVec_(new_ev[c],k)*Ainv_[this->new_c[c]](k,j);
			}
			for(unsigned int i(0);i<this->M_;i++){
				tmp_(i,j) = t_tmp*Ainv_[this->new_c[c]](i,new_r[c])/w[c];
			}
		}
		Ainv_[this->new_c[c]] -= tmp_;
	}
}

template<typename Type>
void SystemFermionic<Type>::swap(){
	MCSystem<Type>::swap();
	new_r[0] = row_(this->new_s[0],this->new_p[0]); 
	new_r[1] = row_(this->new_s[1],this->new_p[1]); 
	new_ev[0] = this->new_c[0]*this->n_ + this->new_s[1];
	new_ev[1] = this->new_c[1]*this->n_ + this->new_s[0];
}

template<typename Type>
void SystemFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	MCSystem<Type>::swap(s0,s1,p0,p1);
	new_r[0] = row_(this->new_s[0],this->new_p[0]); 
	new_r[1] = row_(this->new_s[1],this->new_p[1]); 
	new_ev[0] = this->new_c[0]*this->n_ + this->new_s[1];
	new_ev[1] = this->new_c[1]*this->n_ + this->new_s[0];
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemFermionic<Type>::ratio(){
	if(this->new_c[0] == this->new_c[1]){
		return 1.0;
	} else {
		for(unsigned int c(0);c<2;c++){
			w[c] = 0.0;
			for(unsigned int k(0);k<this->M_;k++){
				w[c] += this->EVec_(new_ev[c],k)*Ainv_[this->new_c[c]](k,new_r[c]);
			}
		}
		/*the minus is correct, it comes from <C|H|C'> because when H is
		 * applied on |C>, the operators are not in the correct color order, so
		 * they need to be exchanged*/
		return -w[0]*w[1];
	}
}
/*}*/
#endif
