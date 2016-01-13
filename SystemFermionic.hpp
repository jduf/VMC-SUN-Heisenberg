#ifndef DEF_SYSTEMFERMIONIC
#define DEF_SYSTEMFERMIONIC

#include "MCSystem.hpp"
#include "Fermionic.hpp"

/*{*//*!Class that makes the connection between MCSystem and Fermionic

	   All MCSystem's pure virtual methods are redefined in this class. It has
	   access to all the necessary informations to sample the configurations of
	   a fermionic system.*//*}*/
template<typename Type>
class SystemFermionic: public MCSystem, public Fermionic<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemFermionic(Fermionic<Type> const& F);
		/*!Constructor that reads from file*/
		SystemFermionic(IOFiles& r);
		/*!Destructor that deletes Ainv and tmp*/
		~SystemFermionic();
		/*{Forbidden*/
		SystemFermionic() = delete;
		SystemFermionic(SystemFermionic<Type>&&) = delete;
		SystemFermionic& operator=(SystemFermionic<Type>) = delete;
		/*}*/

		/*!Set row_ and new_ev_*/
		void swap();
		/*!Set row_ and new_ev_*/
		void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);

		/*{*//*!Computes the ratio of two determinants
			   related to the current and next configuration
			   + when particle of the same color are exchanged a minus sign
			   arises to conserve the Marshall-Peierls sign rule
			   + when two different colors are exchanged, computes the ratio
			   using the determinant lemma
			   *//*}*/
		double ratio();
		/*{*//*!Updates the variables to the new configuration
			   by calling System::update() and then
			   + updates row_
			   + updates the Ainv_ matrices
			   *//*}*/
		void update();

		/*!Returns a copy of this instance of SystemFermionic*/
		std::unique_ptr<MCSystem> clone() const;
		/*!Sets most of the matrices to NULL*/
		void free_memory();
		/*!Writes the current state of the system, the color configuration, the
		 * observables and everything relevant to the simulation*/
		void write(IOFiles& w) const;
		/*!Prints some info related to the system*/
		void print();

	private:
		/*!Authorizes copy only via clone()*/
		SystemFermionic(SystemFermionic<Type> const& SF);

		Matrix<unsigned int> row_;//!< row of the matrix A that is modified
		Matrix<Type>* Ainv_;	  //!< inverse of A
		Matrix<Type>* tmp_;		  //!< temporary matrix used during the update
		Vector<Type>* tmp_v;	  //!< temporary vector used during the update
		Type w_[2];				  //!< det(W)= d = determinant ratios of <GS|a>/<GS|b>; W=(w11,0;0,w22)
		unsigned int new_r_[2];	  //!< selected rows of the A matrices
		unsigned int new_ev_[2];  //!< selected rows of the EVec matrices

		/*!Returns true if the Ainv_ matrices are invertible*/
		bool are_invertible();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
SystemFermionic<Type>::SystemFermionic(Fermionic<Type> const& F):
	System(F),
	MCSystem(F),
	Fermionic<Type>(F),
	row_(n_,m_),
	Ainv_(new Matrix<Type>[N_]),
	tmp_(new Matrix<Type>[N_]),
	tmp_v(new Vector<Type>[N_])
{
	/*!Initialized class variables*/
	for(unsigned int c(0);c<N_;c++){
		Ainv_[c].set(M_(c),M_(c));
		tmp_[c].set(M_(c),M_(c));
		tmp_v[c].set(M_(c));
	}

	/*!Initialized Ainv_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c(0);
	Vector<unsigned int> row(N_,0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				Ainv_[c](row(c),j) = this->EVec_[c](s,j);
			}
			row_(s,p) = row(c);
			row(c)++;
		}
	}

	/*!Make sure that the matrices Ainv_ are invertible by going to a state of
	 * heigh weight*/
	unsigned int l(0);
	unsigned int TRY_MAX(1e5);
	Matrix<Type>* A;
	A = new Matrix<Type>[N_];
	Vector<Type> det_A(N_,1.0);
	Vector<Type> det_Ainv(N_,1.0);
	for(unsigned int c(0);c<N_;c++){
		A[c] = Ainv_[c];
		det_Ainv(c) = Lapack<Type>(Ainv_[c],true,'G').det();
	}
	while( !are_invertible() && ++l<TRY_MAX ){
		swap();
		for(unsigned int j(0);j<M_(new_c_[0]);j++){
			A[new_c_[0]](new_r_[0],j) = this->EVec_[new_c_[0]](new_s_[1],j);
		}
		for(unsigned int j(0);j<M_(new_c_[1]);j++){
			A[new_c_[1]](new_r_[1],j) = this->EVec_[new_c_[1]](new_s_[0],j);
		}
		/*!Compute the ratio of the determinant of the two states. This is the
		 * brute force method but as for now the inverse matrix is unknown, it
		 * is the only solution*/
		Type d(1.0);
		for(unsigned int c(0);c<2;c++){
			det_A(new_c_[c]) = Lapack<Type>(A[new_c_[c]],true,'G').det();
			d *= det_A(new_c_[c])/det_Ainv(new_c_[c]);
		}
		if( my::norm_squared(d)>1 ){
			det_Ainv = det_A;
			/*update the new state*/
			for(unsigned int j(0);j<M_(new_c_[0]);j++){
				Ainv_[new_c_[0]](new_r_[0],j) = this->EVec_[new_c_[0]](new_s_[1],j);
			}
			for(unsigned int j(0);j<M_(new_c_[1]);j++){
				Ainv_[new_c_[1]](new_r_[1],j) = this->EVec_[new_c_[1]](new_s_[0],j);
			}
			s_(new_s_[0],new_p_[0]) = new_c_[1];
			s_(new_s_[1],new_p_[1]) = new_c_[0];
			row_(new_s_[0],new_p_[0]) = new_r_[1];
			row_(new_s_[1],new_p_[1]) = new_r_[0];
		} else {
			/*restore A_ match previous Ainv_ state*/
			for(unsigned int j(0);j<M_(new_c_[0]);j++){
				A[new_c_[0]](new_r_[0],j) = this->EVec_[new_c_[0]](new_s_[0],j);
			}
			for(unsigned int j(0);j<M_(new_c_[1]);j++){
				A[new_c_[1]](new_r_[1],j) = this->EVec_[new_c_[1]](new_s_[1],j);
			}
		}
	}

	/*!Proceed to the inversion if possible*/
	if(l<TRY_MAX){
		status_--;
		for(unsigned int c(0); c<N_; c++){
			Lapack<Type> inv(Ainv_[c],false,'G');
			inv.inv();
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no initial state found after "<<TRY_MAX<<" attempts"<<std::endl; }

	delete[] A;
}

template<typename Type>
SystemFermionic<Type>::SystemFermionic(SystemFermionic<Type> const& SF):
	System(SF),
	MCSystem(SF),
	Fermionic<Type>(SF),
	row_(SF.row_),
	Ainv_(new Matrix<Type>[N_]),
	tmp_(new Matrix<Type>[N_]),
	tmp_v(new Vector<Type>[N_])
{
	/*!Initialized class variables*/
	for(unsigned int c(0);c<N_;c++){ Ainv_[c].set(M_(c),M_(c)); }

	/*!Initialized Ainv_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c(0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				Ainv_[c](row_(s,p),j) = this->EVec_[c](s,j);
			}
		}
	}
	if(are_invertible()){
		for(unsigned int c(0);c<N_;c++){
			Lapack<Type>(Ainv_[c],false,'G').inv();
			tmp_[c].set(M_(c),M_(c));
			tmp_v[c].set(M_(c));
		}
	} else {
		this->status_++;
		std::cerr<<__PRETTY_FUNCTION__<<" the A matrices are not invertible anymore"<<std::endl;
	}
}

template<typename Type>
SystemFermionic<Type>::SystemFermionic(IOFiles& r):
	System(r),
	MCSystem(r),
	Fermionic<Type>(r),
	row_(r),
	Ainv_(new Matrix<Type>[N_]),
	tmp_(new Matrix<Type>[N_]),
	tmp_v(new Vector<Type>[N_])
{
	/*!Initialized class variables*/
	for(unsigned int c(0);c<N_;c++){ Ainv_[c].set(M_(c),M_(c)); }

	/*!Initialized Ainv_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c(0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				Ainv_[c](row_(s,p),j) = this->EVec_[c](s,j);
			}
		}
	}
	if(are_invertible()){
		for(unsigned int c(0);c<N_;c++){
			tmp_[c].set(M_(c),M_(c));
			tmp_v[c].set(M_(c));
			Lapack<Type>(Ainv_[c],false,'G').inv();
		}
	} else {
		this->status_++;
		std::cerr<<__PRETTY_FUNCTION__<<" the A matrices are not invertible anymore "<<std::endl;
	}
}

template<typename Type>
SystemFermionic<Type>::~SystemFermionic(){
	delete[] Ainv_;
	delete[] tmp_;
	delete[] tmp_v;
}

template<typename Type>
std::unique_ptr<MCSystem> SystemFermionic<Type>::clone() const {
	return std::unique_ptr<SystemFermionic<Type> >(new SystemFermionic<Type>(*this));
}
/*}*/

/*void methods*/
/*{*/
template<typename Type>
void SystemFermionic<Type>::swap(){
	MCSystem::swap();
	new_r_[0] = row_(new_s_[0],new_p_[0]);
	new_r_[1] = row_(new_s_[1],new_p_[1]);
	new_ev_[0] = new_s_[1];
	new_ev_[1] = new_s_[0];
}

template<typename Type>
void SystemFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	MCSystem::swap(s0,s1,p0,p1);
	new_r_[0] = row_(new_s_[0],new_p_[0]);
	new_r_[1] = row_(new_s_[1],new_p_[1]);
	new_ev_[0] = new_s_[1];
	new_ev_[1] = new_s_[0];
}

template<typename Type>
void SystemFermionic<Type>::update(){
	MCSystem::update();
	row_(new_s_[0],new_p_[0]) = new_r_[1];
	row_(new_s_[1],new_p_[1]) = new_r_[0];

	//Type t_tmp;
	//unsigned int c_tmp;
	//for(unsigned int c(0);c<2;c++){
	//c_tmp = new_c_[c];
	//for(unsigned int j(0);j<M_(c_tmp);j++){
	//if(new_r_[c] == j){ t_tmp = -1.0; }
	//else { t_tmp = 0.0; }
	//for(unsigned int k(0);k<M_(c_tmp);k++){
	//t_tmp += this->EVec_[c_tmp](new_ev_[c],k)*Ainv_[c_tmp](k,j);
	//}
	//t_tmp /= w_[c];
	//for(unsigned int i(0);i<M_(c_tmp);i++){
	//tmp_[c_tmp](i,j) = t_tmp*Ainv_[c_tmp](i,new_r_[c]);
	//}
	//}
	//Ainv_[c_tmp] -= tmp_[c_tmp];
	//}

	/*remove tmp_*/
	Type Ainvlk;
	unsigned int c;
	unsigned int M;
	unsigned int r;
	for(unsigned int i(0);i<2;i++){
		c = new_c_[i];
		r = new_r_[i];
		M = M_(c);
		/*!compute u.u^T.Ã.A^(-1) = ((A^(-1))^T.Ã^T.u.u^T)^T*/
		BLAS::gemv('T',M,M,Ainv_[c].ptr(),this->EVec_[c].ptr()+new_ev_[i],this->EVec_[c].row(),tmp_v[c].ptr());
		tmp_v[c](r) -= 1.0;
		for(unsigned int j(0);j<M;j++){
			/*need to save this temporary value because Ainv_ is overwritten*/
			Ainvlk = Ainv_[c](j,r)/w_[i];
			for(unsigned int k(0);k<M;k++){ Ainv_[c](j,k) -= Ainvlk*tmp_v[c](k); }
		}
	}
}

template<typename Type>
void SystemFermionic<Type>::write(IOFiles& w) const {
	System::write(w);
	MCSystem::write(w);
	w<<this->same_wf_;
	if(this->same_wf_){ w<<this->EVec_[0]; }
	else { for(unsigned int c(0);c<N_;c++){ w<<this->EVec_[c]; } }
	w<<row_;
}

template<typename Type>
void SystemFermionic<Type>::free_memory(){
	std::cout<<"free"<<std::endl;
	Ainv_[0].set();
	tmp_[0].set();
	tmp_v[0].set();
	for(unsigned int c(1);c<N_;c++){
		if(this->same_wf_){ this->EVec_[c].set(); }
		Ainv_[c].set();
		tmp_[c].set();
		tmp_v[c].set();
	}

	//delete[] Ainv_;
	//Ainv_ = NULL;
	//delete[] tmp_;
	//tmp_ = NULL;
	//delete[] tmp_v;
	//tmp_v = NULL;
	//if(this->same_wf_){
	//for(unsigned int c(1);c<N_;c++){ this->EVec_[c].set(); }
	//}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
double SystemFermionic<Type>::ratio(){
	if(new_c_[0] == new_c_[1]){
		/*!there is no minus sign because if the same color is inverted, the
		 * matrices will be identical up to the inversion of two columns, this
		 * minus sign is then cancelled by the reordering of the operators */
		return 1.0;
	} else {
		//unsigned int c;
		for(unsigned int i(0);i<2;i++){
			w_[i] = BLAS::dot(M_(new_c_[i]),this->EVec_[new_c_[i]].ptr(),true,this->EVec_[new_c_[i]].row(),new_ev_[i],Ainv_[new_c_[i]].ptr(),false,M_(new_c_[i]),new_r_[i]);
			//c= new_c_[i];
			//w_[i] = 0.0;
			//for(unsigned int k(0);k<M_(c);k++){
			//w_[i] += this->EVec_[c](new_ev_[i],k)*Ainv_[c](k,new_r_[i]);
			//}
		}
		/*!the minus sign is correct, it comes from <C|H|C'> because when H is
		 * applied on |C>, the operators are not in the correct color order,
		 * so they need to be exchanged*/
		return -my::real(w_[0]*w_[1]);
	}
}

template<typename Type>
void SystemFermionic<Type>::print(){
	Matrix<Type>* A(new Matrix<Type>[N_]);
	for(unsigned int c(0);c<N_;c++){ A[c].set(M_(c),M_(c)); }
	unsigned int c(0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				A[c](row_(s,p),j) = this->EVec_[c](s,j);
			}
		}
	}
	for(unsigned int c(0);c<N_;c++){ std::cout<<(A[c]*Ainv_[c]).chop()<<std::endl<<std::endl; }
	delete[] A;
}
/*}*/

/*private method*/
/*{*/
template<typename Type>
bool SystemFermionic<Type>::are_invertible(){
	Vector<int> ipiv;
	double rcn(0.0);
	for(unsigned int c(0);c<N_;c++){
		Lapack<Type> inv(Ainv_[c],true,'G');
		ipiv = inv.is_singular(rcn);
		if(!ipiv.ptr()){ return false; }
	}
	return true;
}
/*}*/
#endif
