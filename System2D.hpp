#ifndef DEF_SYSTEM2D
#define DEF_SYSTEM2D

#include "Gnuplot.hpp"
#include "Lapack.hpp"
#include "Combination.hpp"
#include "Rand.hpp"

template<typename Type>
class System2D: public virtual System{
	public:
		/*!Constructor*/
		System2D(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, unsigned int const& sel0, unsigned int const& sel1);
		/*!Destructor*/
		virtual ~System2D();

	protected:
		Matrix<Type> H_;			//!< matrix used to get the band structure
		Matrix<Type> Tx_;			//!< translation operator along x-axis
		Matrix<Type> Ty_;			//!< translation operator along y-axis
		Matrix<std::complex<double> > evec_;//!< translation operator along y-axis
		Vector<double> e_;			//!< eigenvalue of the Hamiltonian T_
		Vector<double> px_;			//!< eigenvalue of the translation along x
		Vector<double> py_;			//!< eigenvalue of the translation along y
		unsigned int const Lx_;		//!< number of unit cell along the x-axis
		unsigned int const Ly_;		//!< number of unit cell along the y-axis
		unsigned int const spuc_;	//!< site per unit cell
		unsigned int sel_[2];
		Vector<unsigned int>* select_;//!< eigenvalue of the Hamiltonian T_

		void compute_TxTy();
		void compute_band_structure();
		void plot_band_structure();
		void select_eigenvectors();

	private:
		/*!Forbids copy*/
		System2D(System2D const& bs);
		/*!Forbids assignment*/
		System2D& operator=(System2D bs);

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);
};
	
/*{constructors*/
template<typename Type>
System2D<Type>::System2D(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, unsigned int const& sel0, unsigned int const& sel1):
	Lx_(sqrt(Lx*n_/(Ly*spuc))),
	Ly_(sqrt(Ly*n_/(Lx*spuc))),
	spuc_(spuc),
	select_(new Vector<unsigned int>[N_])
{
	sel_[0]= sel0;
	sel_[1]= sel1;
}

template<typename Type>
System2D<Type>::~System2D(){
	if(select_){ delete[] select_;}
}
/*}*/

/*{protected methods*/
template<typename Type>
void System2D<Type>::compute_TxTy(){
	Tx_.set(n_,n_,0);
	Ty_.set(n_,n_,0);
	unsigned int tmp;
	double t(1);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			tmp = spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				Tx_(tmp+k, tmp+k+spuc_) = t;
			}
		}
		tmp = spuc_*((Lx_-1) + j*Lx_);
		for(unsigned int k(0);k<spuc_;k++){
			Tx_(tmp+k,spuc_*j*Lx_ + k) = bc_*t;
		}
	}
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_-1;j++){
			tmp = spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				Ty_(tmp+k, tmp+spuc_*Lx_+k) = t;
			}
		}
		tmp = spuc_*(i + (Ly_-1)*Lx_);
		for(unsigned int k(0);k<spuc_;k++){
			Ty_(tmp+k, spuc_*i+k) = bc_*t;
		}
	}
}

template<typename Type>
void System2D<Type>::compute_band_structure(){
	e_.set(n_);
	Rand rnd(1e4);
	Matrix<Type> M;
	bool degenerate;
	unsigned int iter(0);
	do{
		/*may be optimized : */
		M  = H_;
		M += Tx_*Type(rnd.get()*rnd.get(1e6));
		M += Ty_*Type(rnd.get()*rnd.get(1e6));

		Vector<std::complex<double> > eval;
		Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);
		for(unsigned int i(0);i<n_;i++){
			e_(i) = projection(H_,i).real();
		}

		degenerate = false;
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(i+1);j<n_;j++){
				if(are_equal(eval(i),eval(j),1e-10)){
					degenerate = true;
				}
			}
		}
	} while ( degenerate && iter++<1e2);
	if( iter>1e2 ){ std::cerr<<"void System2D<Type>::compute_band_structure() : degenerate"<<std::endl; }
	else {
		Vector<unsigned int> index;
		Matrix<std::complex<double> > evec_tmp(evec_);
		e_.sort(std::less_equal<double>(),index);
		px_.set(n_);
		py_.set(n_);
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<n_;j++){
				std::swap(evec_(i,j),evec_tmp(i,index(j)));
			}
		}
		for(unsigned int i(0);i<n_;i++){
			px_(i) = log(projection(Tx_,i)).imag();
			py_(i) = log(projection(Ty_,i)).imag();
		}
	}
}

template<typename Type>
void System2D<Type>::plot_band_structure(){
	IOFiles spectrum("spectrum.dat",true);
	for(unsigned int i(0);i<n_;i++){
		spectrum<<px_(i)<<" "<<py_(i)<<" "<<e_(i)<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="splot 'spectrum.dat' u 1:2:3";
	gp.save_file();
}

template<typename Type>
void System2D<Type>::select_eigenvectors(){
	//std::cout<<e_<<std::endl;
	unsigned int iter(0);
	double Px(0.0);
	double Py(0.0);
	for(unsigned int c(0);c<N_;c++){
		unsigned int a(M_(c)-1);
		unsigned int b(M_(c)-1);
		//std::cout<<a<<" "<<b<<std::endl;
		do{b++;} while (b+1<n_ && are_equal(e_(b),e_(b-1)));
		//std::cout<<a<<" "<<b<<std::endl;
		if(b!=M_(c)){ while(a>0 && are_equal(e_(a-1),e_(a))){a--;} }
		//std::cout<<a<<" "<<b<<std::endl;
		//std::cout<<e_.range(a,b)<<std::endl;
		Vector<unsigned int> cnk;
		Combination cbn;
		cbn.set(M_(c)-a,b-a,cnk);
		select_[c].set(M_(c));
		for(unsigned int i(0);i<a;i++){ select_[c](i) = i; }
		//do{ 
		//for(unsigned int i(0);i+a<M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		//double Px(0.0);
		//double Py(0.0);
		//double  e(0.0);
		//for(unsigned int i(0);i<M_(c);i++){
		////Px+= (are_equal(std::abs(px_(a+comb(j))),M_PI,1e-8)?0:px_(a+comb(j)));
		////Py+= (are_equal(std::abs(py_(a+comb(j))),M_PI,1e-8)?0:py_(a+comb(j)));
		//e += e_(select_[c](i));
		//}
		//std::cout<<Px<<" "<<Py<<" "<<e<<std::endl;
		//} while (cbn.next());
		while(iter++<sel_[c] && cbn.next());
		for(unsigned int i(0);i+a<M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		double  e(0.0);
		for(unsigned int i(0);i<M_(c);i++){
			Px+= are_equal(std::abs(px_(select_[c](i))),M_PI,1e-12)?0:px_(select_[c](i));
			Py+= are_equal(std::abs(py_(select_[c](i))),M_PI,1e-12)?0:py_(select_[c](i));
			e +=  e_(select_[c](i));
		}
	}
	if(are_equal(Px,0.) && are_equal(Py,0.)){
		std::cout<<iter<<std::endl;
	}
}
/*}*/

/*{private methods*/
template<typename Type>
std::complex<double> System2D<Type>::projection(Matrix<Type> const& O, unsigned int const& idx){
	std::complex<double> tmp;
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		tmp = 0.0;
		for(unsigned int j(0);j<O.col();j++){
			tmp += O(i,j)*evec_(j,idx);
		}
		out += std::conj(evec_(i,idx))*tmp;
	}
	return out;
}
/*}*/
#endif
