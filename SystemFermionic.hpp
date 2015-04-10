#ifndef DEF_SYSTEMFERMIONIC
#define DEF_SYSTEMFERMIONIC

#include "MCSystem.hpp"
#include "Fermionic.hpp"

/*!Class that contains all the necessary informations to sample the
 * configuration of a fermionic system.*/
template<typename Type>
class SystemFermionic : public Fermionic<Type>, public MCSystem<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemFermionic(Fermionic<Type> const& S);
		/*!Destructor that deletes Ainv and tmp*/
		~SystemFermionic();

		/*!Set row and new_ev_*/
		void swap();
		/*!Set row and new_ev_*/
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

		/*!returns true is Ainv_ matrices are invertibles*/
		bool are_invertible();

		Matrix<unsigned int> row_;//!< row of the matrix A that is modified
		Matrix<Type>* Ainv_;	//!< inverse of A
		Matrix<Type>* tmp_;		//!< temporary matrix used during the update 
		Type w_[2];				//!< det(W)= d = determinant ratios of <GS|a>/<GS|b>; W=(w11,0;0,w22)
		unsigned int new_r_[2];	//!< rows of the Ainv_ matrix that are modified (the rows of the related A matrix are modified)
		unsigned int new_ev_[2];//!< newly selected rows of the EVec matrix
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
SystemFermionic<Type>::SystemFermionic(Fermionic<Type> const& S):
	System(S),
	Fermionic<Type>(S),
	MCSystem<Type>(S),
	row_(this->n_,this->m_),
	Ainv_(new Matrix<Type>[this->N_]),
	tmp_(new Matrix<Type>[this->N_])
{ 
	/*!Initialized class variables*/
	for(unsigned int c(0); c<this->N_; c++){ 
		Ainv_[c].set(this->M_(c),this->M_(c)); 
		tmp_[c].set(this->M_(c),this->M_(c));
	}

	/*!Initialized Ainv_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c_tmp(0);
	Vector<unsigned int> row_tmp(this->N_,0);
	for(unsigned int s(0); s < this->n_; s++){
		for(unsigned int p(0); p < this->m_; p++){
			c_tmp = this->s_(s,p);
			for(unsigned int j(0); j < this->M_(c_tmp); j++){
				Ainv_[c_tmp](row_tmp(c_tmp),j) = this->EVec_[c_tmp](s,j);
			}
			row_(s,p) = row_tmp(c_tmp);
			row_tmp(c_tmp)++;
		}
	}

	/*!Make sure that the matrices Ainv_ are invertible by going to a state of
	 * heigh weight*/
	unsigned int l(0);
	unsigned int TRY_MAX(1e5);
	Matrix<Type>* A;
	A = new Matrix<Type>[this->N_];
	Vector<Type> det_A(this->N_,1.0);
	Vector<Type> det_Ainv(this->N_,1.0);
	for(unsigned int c(0); c<this->N_; c++){
		A[c] = Ainv_[c];
		det_Ainv(c) = Lapack<Type>(Ainv_[c],true,'G').det();
	}
	while( !are_invertible() && ++l<TRY_MAX ){
		swap();
		for(unsigned int j(0);j<this->M_(this->new_c_[0]);j++){
			A[this->new_c_[0]](new_r_[0],j) = this->EVec_[this->new_c_[0]](this->new_s_[1],j);
		}
		for(unsigned int j(0);j<this->M_(this->new_c_[1]);j++){
			A[this->new_c_[1]](new_r_[1],j) = this->EVec_[this->new_c_[1]](this->new_s_[0],j);
		}
		/*!Compute the ratio of the determinant of the two states. This is the
		 * brute force method but as for now the inverse matrix is unknown, it
		 * is the only solution*/
		Type d(1.0);
		for(unsigned int c(0);c<2;c++){ 
			det_A(this->new_c_[c]) = Lapack<Type>(A[this->new_c_[c]],true,'G').det();
			d *= det_A(this->new_c_[c])/det_Ainv(this->new_c_[c]);
		}
		if( my::norm_squared(d)>1 ){
			det_Ainv = det_A;
			/*update the new state*/
			for(unsigned int j(0);j<this->M_(this->new_c_[0]);j++){
				Ainv_[this->new_c_[0]](new_r_[0],j) = this->EVec_[this->new_c_[0]](this->new_s_[1],j);
			}
			for(unsigned int j(0);j<this->M_(this->new_c_[1]);j++){
				Ainv_[this->new_c_[1]](new_r_[1],j) = this->EVec_[this->new_c_[1]](this->new_s_[0],j);
			}
			this->s_(this->new_s_[0],this->new_p_[0]) = this->new_c_[1];
			this->s_(this->new_s_[1],this->new_p_[1]) = this->new_c_[0];
			row_(this->new_s_[0],this->new_p_[0]) = new_r_[1];
			row_(this->new_s_[1],this->new_p_[1]) = new_r_[0];
		} else {
			/*restore the s_ and row_ to match Ainv_ state*/
			for(unsigned int j(0);j<this->M_(this->new_c_[0]);j++){
				A[this->new_c_[0]](new_r_[0],j) = this->EVec_[this->new_c_[0]](this->new_s_[0],j);
			}
			for(unsigned int j(0);j<this->M_(this->new_c_[1]);j++){
				A[this->new_c_[1]](new_r_[1],j) = this->EVec_[this->new_c_[1]](this->new_s_[1],j);
			}
		}
	}

	/*!Proceed to the inversion if possible*/
	if(l<TRY_MAX){
		this->status_--; 
		for(unsigned int c(0); c<this->N_; c++){
			Lapack<Type> inv(Ainv_[c],false,'G');
			inv.inv();
		}
	} else { std::cerr<<"No initial state found after "<<TRY_MAX<<" trials"<<std::endl; }

	delete[] A;
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
	new_r_[0] = row_(this->new_s_[0],this->new_p_[0]); 
	new_r_[1] = row_(this->new_s_[1],this->new_p_[1]); 
	new_ev_[0] = this->new_s_[1];
	new_ev_[1] = this->new_s_[0];
}

template<typename Type>
void SystemFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	MCSystem<Type>::swap(s0,s1,p0,p1);
	new_r_[0] = row_(this->new_s_[0],this->new_p_[0]); 
	new_r_[1] = row_(this->new_s_[1],this->new_p_[1]); 
	new_ev_[0] = this->new_s_[1];
	new_ev_[1] = this->new_s_[0];
}

template<typename Type>
void SystemFermionic<Type>::update(){
	MCSystem<Type>::update();
	row_(this->new_s_[0],this->new_p_[0]) = new_r_[1];
	row_(this->new_s_[1],this->new_p_[1]) = new_r_[0];

	Type t_tmp;
	unsigned int c_tmp;
	for(unsigned int c(0);c<2;c++){
		c_tmp = this->new_c_[c];
		for(unsigned int j(0);j<this->M_(c_tmp);j++){
			if(new_r_[c] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<this->M_(c_tmp);k++){
				t_tmp += this->EVec_[c_tmp](new_ev_[c],k)*Ainv_[c_tmp](k,j);
			}
			t_tmp /= w_[c];
			for(unsigned int i(0);i<this->M_(c_tmp);i++){
				tmp_[c_tmp](i,j) = t_tmp*Ainv_[c_tmp](i,new_r_[c]);
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
	if(this->new_c_[0] == this->new_c_[1]){
		/*!there is no minus sign because if the same color is inverted, the
		 * matrices will be identical up to the invertion of two columns, this
		 * minus sign is then cancelled by the reordering of the operators */
		return 1.0;
	} else {
		unsigned int c_tmp;
		for(unsigned int c(0);c<2;c++){
			c_tmp = this->new_c_[c];
			w_[c] = 0.0;
			for(unsigned int k(0);k<this->M_(c_tmp);k++){
				w_[c] += this->EVec_[c_tmp](new_ev_[c],k)*Ainv_[c_tmp](k,new_r_[c]);
			}
		}
		/*!the minus sign is correct, it comes from <C|H|C'> because when H is
		 * applied on |C>, the operators are not in the correct color order, so
		 * they need to be exchanged*/
		return -w_[0]*w_[1];
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

/*private method*/
/*{*/
template<typename Type>
bool SystemFermionic<Type>::are_invertible(){
	Vector<int> ipiv;
	double rcn(0.0);
	for(unsigned int c(0); c<this->N_; c++){
		Lapack<Type> inv(Ainv_[c],true,'G');
		ipiv = inv.is_singular(rcn);
		if(!ipiv.ptr()){ return false; }
	}
	return true;
}
/*}*/
#endif
