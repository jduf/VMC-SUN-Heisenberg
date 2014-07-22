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
		System2D(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc);
		/*!Destructor*/
		virtual ~System2D(){}

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

		void compute_TxTy();
		void compute_band_structure();
		void plot_band_structure();
		void select();

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
System2D<Type>::System2D(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc):
	Lx_(sqrt(Lx*this->n_/(Ly*spuc))),
	Ly_(sqrt(Ly*this->n_/(Lx*spuc))),
	spuc_(spuc)
{}
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
void System2D<Type>::select(){
	for(unsigned int i(0);i<n_;i++){
		std::cout<<e_(i)<<(i+1==M_(0)?"|":" ");
	}
	std::cout<<std::endl;
	unsigned int a(M_(0)-1);
	unsigned int b(M_(0)-1);
	while(are_equal(e_(b+1),e_(b))){b++;}
	b++;
	std::cout<<M_(0)<<" "<<b<<std::endl;
	if(b!=M_(0)){
		while(are_equal(e_(a-1),e_(a))){a--;}
		std::cout<<M_(0)<<" "<<a<<std::endl;
	}
	std::cout<<e_.range(a,b)<<std::endl;
	Vector<unsigned int> comb;
	Combination cbn(M_(0)-a,b-a,comb);
	do{ 
		std::cout<<comb<<std::endl;
		double Px(0.0);
		double Py(0.0);
		double  e(0.0);
		for(unsigned int j(0);j<a;j++){
			Px+= (are_equal(std::abs(px_(j)),M_PI,1e-8)?0:px_(j));
			Py+= (are_equal(std::abs(py_(j)),M_PI,1e-8)?0:px_(j));
			e += e_(j);
		}
		for(unsigned int j(0);a+j<M_(0);j++){
			Px+= (are_equal(std::abs(px_(a+comb(j))),M_PI,1e-8)?0:px_(a+comb(j)));
			Py+= (are_equal(std::abs(py_(a+comb(j))),M_PI,1e-8)?0:py_(a+comb(j)));
			e += e_(a+comb(j));
		}
		std::cout<<e<<" "<<Px<<" "<<Py<<std::endl;
	} while (cbn.next());
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
