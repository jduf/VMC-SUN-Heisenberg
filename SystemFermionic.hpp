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

		/*Initializes Ainv_, tmp_, tmp_v, EVec_ if this instance was
		 * constructed via clone() or SystemFermionic(IOFiles& r).*/
		void init_after_clone_or_reading();

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

		Type determinant_lemma(unsigned int const& i);
		Type determinant_lemma(unsigned int const& col, unsigned int const& s0, unsigned int const& s1, unsigned int const& r0, unsigned int const& r1);
		void compute_peculiar_observable(Observable& O);
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
		status_=0;
		for(unsigned int c(0); c<N_; c++){
			Lapack<Type> inv(Ainv_[c],false,'G');
			inv.inv();
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no initial configuration found after "<<TRY_MAX<<" attempts"<<std::endl; }

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
{ if(!status_){ status_++; } }

template<typename Type>
SystemFermionic<Type>::SystemFermionic(IOFiles& r):
	System(r),
	MCSystem(r),
	Fermionic<Type>(r),
	row_(r),
	Ainv_(new Matrix<Type>[N_]),
	tmp_(new Matrix<Type>[N_]),
	tmp_v(new Vector<Type>[N_])
{ if(!status_){ status_++; } }

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

template<typename Type>
void SystemFermionic<Type>::init_after_clone_or_reading(){
	/*!Initialized class variables*/
	for(unsigned int c(0);c<N_;c++){ Ainv_[c].set(M_(c),M_(c)); }

	/*!If EVec_ are different for each color, then they have already been defined in Fermionic*/
	if(this->same_wf_){ for(unsigned int c(1);c<N_;c++){ this->EVec_[c] = this->EVec_[0]; } }
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
		status_ = 0;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" the A matrices are not invertible anymore"<<std::endl; }
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

	Type t_tmp;
	unsigned int c;
	for(unsigned int i(0);i<2;i++){
		c = new_c_[i];
		for(unsigned int j(0);j<M_(c);j++){
			if(new_r_[i] == j){ t_tmp = -1.0; }
			else { t_tmp = 0.0; }
			for(unsigned int k(0);k<M_(c);k++){
				t_tmp += this->EVec_[c](new_ev_[i],k)*Ainv_[c](k,j);
			}
			t_tmp /= w_[i];
			for(unsigned int k(0);k<M_(c);k++){
				tmp_[c](k,j) = t_tmp*Ainv_[c](k,new_r_[i]);
			}
		}
		Ainv_[c] -= tmp_[c];
	}

	/*remove tmp_*/
	//Type Ainvlk;
	//unsigned int c;
	//unsigned int M;
	//unsigned int r;
	//for(unsigned int i(0);i<2;i++){
	//c = new_c_[i];
	//r = new_r_[i];
	//M = M_(c);
	///*!compute u.u^T.Ã.A^(-1) = ((A^(-1))^T.Ã^T.u.u^T)^T*/
	//BLAS::gemv('T',M,M,Ainv_[c].ptr(),this->EVec_[c].ptr()+new_ev_[i],this->EVec_[c].row(),tmp_v[c].ptr());
	//tmp_v[c](r) -= 1.0;
	//for(unsigned int j(0);j<M;j++){
	///*need to save this temporary value because Ainv_ is overwritten*/
	//Ainvlk = Ainv_[c](j,r)/w_[i];
	//for(unsigned int k(0);k<M;k++){ Ainv_[c](j,k) -= Ainvlk*tmp_v[c](k); }
	//}
	//}
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
	delete[] Ainv_; Ainv_ = NULL;
	delete[] tmp_;  tmp_  = NULL;
	delete[] tmp_v; tmp_v = NULL;
	if(this->same_wf_){
		for(unsigned int c(1);c<N_;c++){ this->EVec_[c].set(); }
	}
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
		/*!the minus sign is correct, it comes from <C|H|C'> because when H is
		 * applied on |C>, the operators are not in the correct color order,
		 * so they need to be exchanged*/
		return -my::real(determinant_lemma(0)*determinant_lemma(1));
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
Type SystemFermionic<Type>::determinant_lemma(unsigned int const& i){
	unsigned int c(new_c_[i]);
	//w_[i] = BLAS::dot(M_(new_c_[i]),this->EVec_[new_c_[i]].ptr(),true,this->EVec_[new_c_[i]].row(),new_ev_[i],Ainv_[new_c_[i]].ptr(),false,M_(new_c_[i]),new_r_[i]);
	w_[i] = 0.0;
	for(unsigned int k(0);k<M_(c);k++){
		w_[i] += this->EVec_[c](new_ev_[i],k)*Ainv_[c](k,new_r_[i]);
	}
	return w_[i];
}

template<typename Type>
Type SystemFermionic<Type>::determinant_lemma(unsigned int const& col, unsigned int const& s0, unsigned int const& s1, unsigned int const& r0, unsigned int const& r1){
	Type a(0.0);
	Type b(0.0);
	Type c(0.0);
	Type d(0.0);
	unsigned int M(this->M_(col));
	for(unsigned int k(0);k<M;k++){
		a += this->EVec_[col](s1,k)*Ainv_[col](k,r0);
		c += this->EVec_[col](s0,k)*Ainv_[col](k,r0);
		c += this->EVec_[col](s1,k)*Ainv_[col](k,r1);
		d += this->EVec_[col](s0,k)*Ainv_[col](k,r1);
	}
	return a*d-b*c;
}

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

template<typename Type>
void SystemFermionic<Type>::compute_peculiar_observable(Observable& O){
	if(O.get_type()==4){
		Data<double>* H2(&O[0]);
		H2->set_x(0.0);

		Type w;
		unsigned int L(this->obs_[0].nlinks());

		auto colors_involved = [](unsigned int const& ci0, unsigned int const& ci1, unsigned int const& cj0, unsigned int const& cj1){
			if(ci0 == ci1){
				if(ci0 == cj0){
					if(ci0 == cj1){ /*AAAA*/ return 10; }
					else { /*AAAB*/ return 30; }
				} else {
					if(ci0 == cj1){ /*AABA*/ return 31; }
					else if( cj0 == cj1) {  /*AABB*/ return 11; }
					else {  /*AABC*/ return 32; }
				}
			} else if(ci0 == cj0){
				if(ci0 == cj1){ /*ABAA*/ return 20; }
				else if (ci1 == cj1) { /*ABAB*/ return 50; }
				else { /*ABAC*/ return 40; }
			} else if(ci0 == cj1){
				if(ci1 == cj0){ /*ABBA*/ return 51; }
				else {  /*ABCA*/ return 41; }
			} else {
				if(ci1 == cj0){
					if(ci1 == cj1){ /*ABBB*/ return 21; }
					else { /*ABBC*/ return 42; }
				} else {
					if(cj0 == cj1){ /*ABCC*/ return 22; }
					else if ( ci1 == cj1){  /*ABCB*/ return 43; }
					else {  /*ABCD*/ return 60; }
				}
			}
		};

		int* idx[2];
		idx[0] = this->obs_[0].get_links().ptr();
		idx[1] = idx[0]+L;
		unsigned int ci0;
		unsigned int ci1;
		unsigned int cj0;
		unsigned int cj1;
		for(unsigned int i(0);i<L;i++){
			int* jdx[2];
			jdx[0] = this->obs_[0].get_links().ptr();
			jdx[1] = jdx[0]+L;
			for(unsigned int j(0);j<L;j++){
				if(
						i>j && /*!To avoid the double computation*/
						*idx[0] != *jdx[0] &&
						*idx[0] != *jdx[1] &&
						*idx[1] != *jdx[0] &&
						*idx[1] != *jdx[1]
				  ){
					/*!Consider here only pairs of links that do not share any
					 * site, they can therefore be treated independently*/
					for(unsigned int ip0(0);ip0<this->m_;ip0++){
						ci0 = this->s_(*idx[0],ip0);
						for(unsigned int ip1(0);ip1<this->m_;ip1++){
							ci1 = this->s_(*idx[1],ip1);
							swap(*idx[0],*idx[1],ip0,ip1);

							if(!this->is_new_state_forbidden()){
								for(unsigned int jp0(0);jp0<this->m_;jp0++){
									cj0 = this->s_(*jdx[0],jp0);
									for(unsigned int jp1(0);jp1<this->m_;jp1++){
										cj1 = this->s_(*jdx[1],jp1);
										swap(*jdx[0],*jdx[1],jp0,jp1);

										if(!this->is_new_state_forbidden()){
											switch(colors_involved(ci0,ci1,cj0,cj1)){
												case 60:/*ABCD*/
													{
														swap(*idx[0],*idx[1],ip0,ip1);
														w = ratio();
														swap(*jdx[0],*jdx[1],jp0,jp1);
														w*= ratio();
														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 51:/*ABBA*/
													{
														w = determinant_lemma(ci0,*idx[0],*jdx[1],row_(*idx[0],ip0),row_(*jdx[1],jp1));
														w*= determinant_lemma(ci1,*idx[1],*jdx[0],row_(*idx[1],ip1),row_(*jdx[0],jp0));
														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 50:/*ABAB*/
													{
														w = determinant_lemma(ci0,*idx[0],*jdx[0],row_(*idx[0],ip0),row_(*jdx[0],jp0));
														w*= determinant_lemma(ci1,*idx[1],*jdx[1],row_(*idx[1],ip1),row_(*jdx[1],jp1));
														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 43:/*ABCB*/
													{
														swap(*idx[0],*idx[1],ip0,ip1);
														w = determinant_lemma(0);
														swap(*jdx[0],*jdx[1],ip0,ip1);
														w*= determinant_lemma(0);
														w*= determinant_lemma(ci1,*idx[1],*jdx[1],row_(*idx[1],ip1),row_(*jdx[1],jp1));

														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 42:/*ABBC*/
													{
														swap(*idx[0],*idx[1],ip0,ip1);
														w = determinant_lemma(1);
														swap(*jdx[0],*jdx[1],ip0,ip1);
														w*= determinant_lemma(1);
														w*= determinant_lemma(ci1,*idx[1],*jdx[0],row_(*idx[1],ip1),row_(*jdx[0],jp0));

														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 41:/*ABCA*/
													{
														swap(*idx[0],*idx[1],ip0,ip1);
														w = determinant_lemma(0);
														swap(*jdx[0],*jdx[1],ip0,ip1);
														w*= determinant_lemma(1);
														w*= determinant_lemma(ci0,*idx[0],*jdx[1],row_(*idx[0],ip0),row_(*jdx[1],jp1));

														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 40:/*ABAC*/
													{
														swap(*idx[0],*idx[1],ip0,ip1);
														w = determinant_lemma(1);
														swap(*jdx[0],*jdx[1],ip0,ip1);
														w*= determinant_lemma(0);
														w*= determinant_lemma(ci0,*idx[0],*jdx[0],row_(*idx[0],ip0),row_(*jdx[0],jp0));

														H2->add(2*J_(i)*J_(j)*my::real(w));
													}break;
												case 32:/*AABC*/
												case 31:/*AABA*/
												case 30:/*AAAB*/
													{
														swap(*idx[0],*idx[1],ip0,ip1);
														H2->add(2*J_(i)*J_(j)*ratio());
													}break;
												case 22:/*ABCC*/
												case 21:/*ABBB*/
												case 20:/*ABAA*/
													{
														swap(*jdx[0],*jdx[1],jp0,jp1);
														H2->add(2*J_(i)*J_(j)*ratio());
													}break;
												case 11:/*AABB*/
												case 10:/*AAAA*/
													{ H2->add(2*J_(i)*J_(j)); }break;
											}
										}
									}
								}
							}
						}
					}
				} else if (i==j){
					H2->add(m_*m_*J_(i)*J_(j));
				} else {
					/*!Consider here only pairs of links that share only one
					 * site which will always be labelled by b. The operator H2
					 * will act on three sites a,b,c and can be decomposed in a
					 * product of two operators Hab and Hbc. */
					unsigned int a(0);
					unsigned int b(0);
					unsigned int c(0);
					if     ( *idx[0] == *jdx[0] ){ a = *idx[1]; b = *idx[0]; c = *jdx[1]; }
					else if( *idx[0] == *jdx[1] ){ a = *idx[1]; b = *idx[0]; c = *jdx[0]; }
					else if( *idx[1] == *jdx[0] ){ a = *idx[0]; b = *idx[1]; c = *jdx[1]; }
					else if( *idx[1] == *jdx[1] ){ a = *idx[0]; b = *idx[1]; c = *jdx[0]; }

					if(a || b || c){
						/*!Thanks to the definition of a,b,c the common site for
						 * the two permutation is b, therefore for example when
						 * a state is symbolised by ABCD, it corresponds to
						 * abbc and after applying Hab*Hbc it becomes BADC
						 * because on the site B there are two different
						 * particles BC. If the initial state in ABBC, then
						 * Hab*Hbc applies on the same particle on the site B
						 * which implies that after the application of Hbc, it
						 * becomes ACCB and after Hab it becomes CAAB.*/
						for(unsigned int ip0(0);ip0<this->m_;ip0++){
							ci0 = this->s_(a,ip0);
							for(unsigned int ip1(0);ip1<this->m_;ip1++){
								ci1 = this->s_(b,ip1);
								for(unsigned int jp0(0);jp0<this->m_;jp0++){
									cj0 = this->s_(b,jp0);
									for(unsigned int jp1(0);jp1<this->m_;jp1++){
										cj1 = this->s_(c,jp1);
										/*test authorised states*/
										switch(colors_involved(ci0,ci1,cj0,cj1)){
											case 60:/*ABCD*/
												{
													bool allowed(true);
													for(unsigned int k(0);k<this->m_;k++){
														if(this->s_(a,k) == ci1 || this->s_(b,k) == ci0 || this->s_(b,k) == cj1 || this->s_(c,k) == cj0) { allowed = false; k=this->m_; }
													}
													if(allowed){ 
														swap(*idx[0],*idx[1],ip0,ip1);
														w = ratio();/*det(A)det(B)*/
														swap(*jdx[0],*jdx[1],jp0,jp1);
														w*= ratio();/*det(C)det(D)*/
														H2->add(J_(i)*J_(j)*my::real(w));
													}
												}break;
											case 51:/*ABBA*/
												{
													swap(b,c,jp0,jp1);
													if(!this->is_new_state_forbidden()){ H2->add(J_(i)*J_(j)*ratio()); }
												}break;
											case 50:/*ABAB*/
												{
													swap(a,c,ip0,jp1);
													if(!this->is_new_state_forbidden()){ H2->add(J_(i)*J_(j)*ratio()); }
												}break;
											case 43:/*ABCB*/
												{
													bool allowed(true);
													for(unsigned int k(0);k<this->m_;k++){
														if(this->s_(a,k) == ci1 || this->s_(b,k) == ci0 || this->s_(c,k) == cj0) { allowed = false; k=this->m_; }
													}
													if(allowed){ 
														swap(a,b,ip0,ip1);
														w = determinant_lemma(0);
														w*= determinant_lemma(1);
														swap(b,c,jp0,jp1);
														w*= determinant_lemma(1);
														H2->add(J_(i)*J_(j)*my::real(w));
													}
												}break;
											case 42:/*ABBC*/
												{
													bool allowed(true);
													for(unsigned int k(0);k<this->m_;k++){
														if(this->s_(a,k) == cj1 || this->s_(b,k) == ci0 || this->s_(c,k) == ci1) { allowed = false; k=this->m_; }
													}
													if(allowed){
														swap(a,b,ip0,ip1);
														w = determinant_lemma(0);/*det(A)*/
														swap(b,c,ip1,jp1);
														w*= determinant_lemma(1);/*det(B)*/
														swap(a,c,ip0,jp1);
														w*= determinant_lemma(1);/*det(C)*/
														H2->add(J_(i)*J_(j)*my::real(w));
													}
												}break;
											case 41:/*ABCA*/
												{ /*impossible output state*/ }break;
											case 40:/*ABAC*/
												{
													bool allowed(true);
													for(unsigned int k(0);k<this->m_;k++){
														if(this->s_(a,k) == ci1 || this->s_(b,k) == cj1 || this->s_(c,k) == cj0) { allowed = false; k=this->m_; }
													}
													if(allowed){
														swap(a,b,ip0,ip1);
														w = determinant_lemma(0);
														w*= determinant_lemma(1);
														swap(b,c,jp0,jp1);
														w*= determinant_lemma(1);
														H2->add(J_(i)*J_(j)*my::real(w));
													}
												}break;
											case 32:/*AABC*/
												{
													swap(b,c,jp0,jp1);
													if(!this->is_new_state_forbidden()){ H2->add(J_(i)*J_(j)*ratio()); }
												}break;
											case 31:/*AABA*/
												{ /*impossible output state*/ }break;
											case 30:/*AAAB*/
												{
													swap(a,c,ip0,jp1);
													if(!this->is_new_state_forbidden()){ H2->add(J_(i)*J_(j)*ratio()); }
												}break;
											case 22:/*ABCC*/
												{
													swap(a,b,ip0,ip1);
													if(!this->is_new_state_forbidden()){ H2->add(J_(i)*J_(j)*ratio()); }
												}break;
											case 21:/*ABBB*/
												{
													swap(a,b,ip0,ip1);
													if(!this->is_new_state_forbidden()){ H2->add(J_(i)*J_(j)*ratio()); }
												}break;
											case 20:/*ABAA*/
												{ /*impossible output state*/ }break;
											case 11:/*AABB*/
											case 10:/*AAAA*/
												{ H2->add(J_(i)*J_(j)); }break;
										}
									}
								}
							}
						}
					}
				}
				jdx[0]++;
				jdx[1]++;
			}
			idx[0]++;
			idx[1]++;
		}
		H2->divide(this->n_*this->n_);
	}
}
/*}*/
#endif
