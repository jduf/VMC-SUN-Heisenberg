#ifndef DEF_SYSTEMBIFERMIONIC
#define DEF_SYSTEMBIFERMIONIC

#include "MCSystem.hpp"
#include "BiFermionic.hpp"

/*!Class that contains all the necessary informations to sample the
 * configuration of a fermionic system.*/
template<typename Type>
class SystemBiFermionic: public MCSystem, public BiFermionic<Type>{
	public:
		/*!Constructor that creates an initial state*/
		SystemBiFermionic(Fermionic<Type> const& F0, Fermionic<Type> const& F1);
		/*!Constructor that reads from file*/
		SystemBiFermionic(IOFiles& r);
		/*!Destructor that deletes Ainv and tmp*/
		~SystemBiFermionic();
		/*{Forbidden*/
		SystemBiFermionic() = delete;
		SystemBiFermionic(SystemBiFermionic<Type>&&) = delete;
		SystemBiFermionic& operator=(SystemBiFermionic<Type>) = delete;
		/*}*/

		/*!Set row_ and new_ev_*/
		void swap();
		/*!Set row_ and new_ev_*/
		void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);

		/*{Description
		 *!Computes the ratio of the two determinants related to the current
		 * and next configuration
		 * - when particle of the same color are exchanged a minus sign arises
		 *   to conserve the Marshall-Peierls sign rule
		 * - when two different colors are exchanged, computes the ratio using
		 *   the determinant lemma
		 }*/
		double ratio();
		/*{Description*/
		/*!Calls System::update() and then
		 * - updates row_
		 * - updates the A_ matrices */
		/*}*/
		void update();
		/*!Measures the system for the new step*/
		void measure_new_step();

		/*!Returns a copy of this instance of SystemBiFermionic*/
		std::unique_ptr<MCSystem> clone() const;
		/*!Sets most of the matrices to NULL*/
		void free_memory();
		/*!Writes the curent state of the system, the color configuration, the
		 * observables and everything relevent to the simulation*/
		void write(IOFiles& w) const;

	private:
		/*!Authorizes copy only via clone()*/
		SystemBiFermionic(SystemBiFermionic<Type> const& SBF);

		bool ratio_for_measure_;  //!< allows the correct matrices' selection to compute the ratio
		Matrix<unsigned int> row_;//!< row of the matrix A that is modified
		Matrix<Type>* A_[2];	  //!< A matrices (one for each wavefunction and colors)
		Matrix<Type> det_;		  //!< determinant of the A matrices
		unsigned int new_r_[2];	  //!< selected rows of the A matrices
		unsigned int new_ev_[2];  //!< selected rows of the EVec matrices
		unsigned int oi_; 		  //!< overlap index : index of the observable related to <0|1>
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
SystemBiFermionic<Type>::SystemBiFermionic(Fermionic<Type> const& F0, Fermionic<Type> const& F1):
	System(F0),
	MCSystem(F0),
	BiFermionic<Type>(F0,F1),
	ratio_for_measure_(true),
	row_(n_,m_),
	A_{new Matrix<Type>[N_],new Matrix<Type>[N_]},
	det_(2,N_),
	oi_(obs_.size())
{
	/*!Initialized class variables*/
	obs_.push_back(Observable("Overlap",10,2,0));
	for(unsigned int c(0);c<N_;c++){
		A_[0][c].set(M_(c),M_(c));
		A_[1][c].set(M_(c),M_(c));
	}

	/*!Initialized A_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c(0);
	Vector<unsigned int> row(N_,0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				A_[0][c](row(c),j) = this->EVec_[0][c](s,j);
				A_[1][c](row(c),j) = this->EVec_[1][c](s,j);
			}
			row_(s,p) = row(c);
			row(c)++;
		}
	}

	for(unsigned int c(0);c<N_;c++){
		det_(0,c) = Lapack<Type>(A_[0][c],true,'G').det();
		det_(1,c) = Lapack<Type>(A_[1][c],true,'G').det();
	}

	if(status_){ status_=0; }
}

template<typename Type>
SystemBiFermionic<Type>::SystemBiFermionic(SystemBiFermionic<Type> const& SBF):
	System(SBF),
	MCSystem(SBF),
	BiFermionic<Type>(SBF),
	ratio_for_measure_(true),
	row_(SBF.row_),
	A_{new Matrix<Type>[N_],new Matrix<Type>[N_]},
	det_(SBF.det_),
	oi_(0)
{
	/*!Initialized class variables*/
	for(unsigned int i(0);i<obs_.size();i++){
		if(obs_[i].get_type() == 10){ oi_ = i; i = obs_.size(); }
	}
	for(unsigned int c(0);c<N_;c++){
		A_[0][c].set(M_(c),M_(c));
		A_[1][c].set(M_(c),M_(c));
	}

	/*!Initialized A_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c(0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				A_[0][c](row_(s,p),j) = this->EVec_[0][c](s,j);
				A_[1][c](row_(s,p),j) = this->EVec_[1][c](s,j);
			}
		}
	}
}

template<typename Type>
SystemBiFermionic<Type>::SystemBiFermionic(IOFiles& r):
	System(r),
	MCSystem(r),
	BiFermionic<Type>(r),
	ratio_for_measure_(true),
	row_(r),
	A_{new Matrix<Type>[N_],new Matrix<Type>[N_]},
	det_(r)
{
	/*!Initialized class variables*/
	for(unsigned int i(0);i<obs_.size();i++){
		if(obs_[i].get_type() == 10){ oi_ = i; i = obs_.size(); }
	}
	for(unsigned int c(0);c<N_;c++){
		A_[0][c].set(M_(c),M_(c));
		A_[1][c].set(M_(c),M_(c));
	}

	/*!Initialized A_ and row_ with the correct eigenvectors according to s_*/
	unsigned int c(0);
	for(unsigned int s(0);s<n_;s++){
		for(unsigned int p(0);p<m_;p++){
			c = s_(s,p);
			for(unsigned int j(0);j<M_(c);j++){
				A_[0][c](row_(s,p),j) = this->EVec_[0][c](s,j);
				A_[1][c](row_(s,p),j) = this->EVec_[1][c](s,j);
			}
		}
	}
}

template<typename Type>
SystemBiFermionic<Type>::~SystemBiFermionic(){
	delete[] A_[0];
	delete[] A_[1];
}

template<typename Type>
std::unique_ptr<MCSystem> SystemBiFermionic<Type>::clone() const {
	return std::unique_ptr<SystemBiFermionic<Type> >(new SystemBiFermionic<Type>(*this));
}
/*}*/

/*void methods*/
/*{*/
template<typename Type>
void SystemBiFermionic<Type>::swap(){
	MCSystem::swap();
	new_r_[0] = row_(new_s_[0],new_p_[0]);
	new_r_[1] = row_(new_s_[1],new_p_[1]);
	new_ev_[0] = new_s_[1];
	new_ev_[1] = new_s_[0];
}

template<typename Type>
void SystemBiFermionic<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	MCSystem::swap(s0,s1,p0,p1);
	new_r_[0] = row_(new_s_[0],new_p_[0]);
	new_r_[1] = row_(new_s_[1],new_p_[1]);
	new_ev_[0] = new_s_[1];
	new_ev_[1] = new_s_[0];
}

template<typename Type>
void SystemBiFermionic<Type>::update(){
	MCSystem::update();
	row_(new_s_[0],new_p_[0]) = new_r_[1];
	row_(new_s_[1],new_p_[1]) = new_r_[0];

	unsigned int c;
	unsigned int r;
	unsigned int ev;
	for(unsigned int i(0);i<2;i++){
		c = new_c_[i];
		r = new_r_[i];
		ev= new_ev_[i];
		for(unsigned int j(0);j<M_(i);j++){
			A_[0][c](r,j) = this->EVec_[0][c](ev,j);
			A_[1][c](r,j) = this->EVec_[1][c](ev,j);
		}
		det_(0,c) = Lapack<Type>(A_[0][c],true,'G').det();
		det_(1,c) = Lapack<Type>(A_[1][c],true,'G').det();
	}
}

template<typename Type>
void SystemBiFermionic<Type>::measure_new_step(){
	ratio_for_measure_ = true;
	MCSystem::measure_new_step();

	Type r(1.0);
	for(unsigned int c(0);c<N_;c++){
		r *= det_(1,c)/det_(0,c);
	}
	obs_[oi_][0].set_x(my::real(r));
	obs_[oi_][1].set_x(my::imag(r));

	ratio_for_measure_ = false;
}

template<typename Type>
void SystemBiFermionic<Type>::write(IOFiles& w) const {
	System::write(w);
	MCSystem::write(w);
	w<<this->same_wf_[0]<<this->same_wf_[1];
	for(unsigned int i(0);i<2;i++){
		if(this->same_wf_[i]){ w<<this->EVec_[i][0]; }
		else { for(unsigned int c(0);c<N_;c++){ w<<this->EVec_[i][c]; } }
	}
	w<<row_<<det_;
}

template<typename Type>
void SystemBiFermionic<Type>::free_memory(){
	std::cout<<"free"<<std::endl;
	A_[0][0].set();
	A_[1][0].set();
	for(unsigned int c(1);c<N_;c++){
		if(this->same_wf_[0]){ this->EVec_[0][c].set(); }
		if(this->same_wf_[1]){ this->EVec_[1][c].set(); }
		A_[0][c].set();
		A_[1][c].set();
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
double SystemBiFermionic<Type>::ratio(){
	/*!the ratio computed to find the new configuration only involve one
	 * wavefunction, hence m=0. for the measurements, the ratio invoves
	 * two wavefunctions, hence m=1 */
	unsigned int m(ratio_for_measure_?1:0);
	Type r(1.0);
	if(new_c_[0] == new_c_[1]){
		if(m){
			for(unsigned int c(0);c<N_;c++){ r *= det_(m,c)/det_(0,c); }
			/*!there is no minus sign because if the same color is inverted,
			 * the matrices will be identical up to the inversion of two
			 * columns, this minus sign is then cancelled by the reordering of
			 * the operators */
			return my::real(r);
		} else { return 1; }
	} else {
		Matrix<Type> Atmp;
		for(unsigned int c(0);c<N_;c++){
			if(c == new_c_[0] || c == new_c_[1]){
				unsigned int i(c==new_c_[0]?0:1);
				Atmp = A_[m][c];
				for(unsigned int j(0);j<M_(c);j++){
					Atmp(new_r_[i],j) = this->EVec_[m][c](new_ev_[i],j);
				}
				r *= Lapack<Type>(Atmp,false,'G').det()/det_(0,c);
			} else if(m){ r *= det_(m,c)/det_(0,c); }
		}
		/*!the minus sign is correct, it comes from <C|H|C'> because when H is
		 * applied on |C>, the operators are not in the correct color order, so
		 * they need to be exchanged*/
		return -my::real(r);
	}
}
/*}*/
#endif
