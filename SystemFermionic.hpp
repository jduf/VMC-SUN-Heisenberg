#ifndef DEF_SYSTEMFERMIONIC
#define DEF_SYSTEMFERMIONIC

#include "MCSystem.hpp"

/*!Class that contains all the necessary informations to sample the configuration of a fermionic system.*/
template<typename Type>
class SystemFermionic : public Fermionic<Type>, public MCSystem<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemFermionic(Fermionic<Type> const& S, Rand& seed);
		/*!Destructor that deletes Ainv and tmp*/
		~SystemFermionic();

		/*!Set row and new_ev*/
		void swap();
		/*!Set row and new_ev*/
		void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);

		/*{Description
		*!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 * - when particle of the same color are exchanged a minus sign arises 
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio using
		 *   the determinant lemma
		}*/
		Type ratio();
		/*{Description*/
		/*!Calls System::update() and then
		 * - updates row_
		 * - updates the Ainv_ matrices */
		/*}*/
		void update();
		/*!Prints some info related to the system*/
		void print();

	private:
		/*!Forbids copy*/
		SystemFermionic(SystemFermionic const& S);
		/*!Forbids assignment*/
		SystemFermionic& operator=(SystemFermionic const& S);

		Matrix<unsigned int> row_;
		Matrix<Type> *Ainv_;	//!< inverse of A
		Matrix<Type>* tmp_;		//!< temporary matrix used during the update 
		Type w[2];				//!< det(W)= d = determinant ratios of <GS|a>/<GS|b>; W=(w11,0,0,w22)
		unsigned int new_r[2];	//!< rows of the Ainv_ matrix that are modified (the rows of the related A matrix are modified)
		unsigned int new_ev[2]; //!< newly selected rows of the EVec matrix
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
SystemFermionic<Type>::SystemFermionic(Fermionic<Type> const& S, Rand& seed):
	System(S),
	Fermionic<Type>(S),
	MCSystem<Type>(S,seed)
{ 
	Ainv_ = new Matrix<Type>[this->N_];
	tmp_ = new Matrix<Type>[this->N_];
	for(unsigned int c(0); c<this->N_; c++){ 
		Ainv_[c].set(this->M_(c),this->M_(c)); 
		tmp_[c].set(this->M_(c),this->M_(c));
	}
	row_.set(this->n_,this->m_);

	unsigned int c(0);
	Vector<unsigned int> M_tmp(this->M_);
	for(unsigned int p(0); p<this->m_; p++){
		for(unsigned int s(0); s<this->n_; s++){
			this->s_(s,p) = c;
			M_tmp(c) -= 1;
			if(!M_tmp(c)){ c++; }
		}
	}

	Vector<int> ipiv;
	unsigned int TRY_MAX(100);
	unsigned int l(0);
	double rcn(0.0);
	do {
		for(unsigned int i(0);i<this->N_*(this->n_*this->m_)*(this->n_*this->m_);i++){
			MCSystem<Type>::swap();
			this->s_(this->new_s[0],this->new_p[0]) = this->new_c[1];
			this->s_(this->new_s[1],this->new_p[1]) = this->new_c[0];
		}
		Vector<unsigned int> row_tmp(this->N_,0);
		for(unsigned int s(0); s < this->n_; s++){
			for(unsigned int p(0); p < this->m_; p++){
				c = this->s_(s,p);
				for(unsigned int j(0); j < this->M_(c); j++){
					Ainv_[c](row_tmp(c),j) = this->EVec_[c](s,j);
				}
				row_(s,p) = row_tmp(c);
				row_tmp(c)++;
			}
		}
		for(unsigned int i(0); i<this->N_; i++){
			Lapack<Type> inv(Ainv_[i],false,'G');
			ipiv = inv.is_singular(rcn);
			if(ipiv.ptr()){ inv.inv(ipiv); }
			else { i = this->N_; }
		}
	} while (!ipiv.ptr() && ++l<TRY_MAX);
	if(l!=TRY_MAX){ this->status_--; }
	else { std::cerr<<"No initial state found after "<<TRY_MAX<<"trials"<<std::endl; }
}

template<typename Type>
SystemFermionic<Type>::~SystemFermionic(){
	delete[] Ainv_;
	delete[] tmp_;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void SystemFermionic<Type>::swap(){
	MCSystem<Type>::swap();
	new_r[0] = row_(this->new_s[0],this->new_p[0]); 
	new_r[1] = row_(this->new_s[1],this->new_p[1]); 
	new_ev[0] = this->new_s[1];
	new_ev[1] = this->new_s[0];
}

template<typename Type>
void SystemFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	MCSystem<Type>::swap(s0,s1,p0,p1);
	new_r[0] = row_(this->new_s[0],this->new_p[0]); 
	new_r[1] = row_(this->new_s[1],this->new_p[1]); 
	new_ev[0] = this->new_s[1];
	new_ev[1] = this->new_s[0];
}

template<typename Type>
void SystemFermionic<Type>::update(){
	MCSystem<Type>::update();
	row_(this->new_s[0],this->new_p[0]) = new_r[1];
	row_(this->new_s[1],this->new_p[1]) = new_r[0];

	Type t_tmp;
	unsigned int c_tmp;
	for(unsigned int c(0);c<2;c++){
		c_tmp = this->new_c[c];
		for(unsigned int j(0);j<this->M_(c_tmp);j++){
			if(new_r[c] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<this->M_(c_tmp);k++){
				t_tmp += this->EVec_[c_tmp](new_ev[c],k)*Ainv_[c_tmp](k,j);
			}
			t_tmp /= w[c];
			for(unsigned int i(0);i<this->M_(c_tmp);i++){
				tmp_[c_tmp](i,j) = t_tmp*Ainv_[c_tmp](i,new_r[c]);
			}
		}
		Ainv_[c_tmp] -= tmp_[c_tmp];
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type> 
Type SystemFermionic<Type>::ratio(){
	if(this->new_c[0] == this->new_c[1]){
		return 1.0;
	} else {
		unsigned int c_tmp;
		for(unsigned int c(0);c<2;c++){
			c_tmp = this->new_c[c];
			w[c] = 0.0;
			for(unsigned int k(0);k<this->M_(c_tmp);k++){
				w[c] += this->EVec_[c_tmp](new_ev[c],k)*Ainv_[c_tmp](k,new_r[c]);
			}
		}
		/*!the minus is correct, it comes from <C|H|C'> because when H is
		 * applied on |C>, the operators are not in the correct color order, so
		 * they need to be exchanged*/
		return -w[0]*w[1];
	}
}

template<typename Type>
void SystemFermionic<Type>::print(){
	Matrix<Type>* A(new Matrix<Type>[this->N_]);
	for(unsigned int c(0);c<this->N_;c++){ A[c].set(this->M_(c),this->M_(c)); } 
	unsigned int c(0);
	for(unsigned int s(0);s<this->n_;s++){
		for(unsigned int p(0);p<this->m_;p++){
			c = this->s_(s,p);
			for(unsigned int j(0);j<this->M_(c);j++){
				A[c](row_(s,p),j) = this->EVec_[c](s,j);
			}
		}
	}
	for(unsigned int c(0);c<this->N_;c++){ std::cout<<(A[c]*Ainv_[c]).chop()<<std::endl<<std::endl; }
	delete[] A;
}
/*}*/
#endif
