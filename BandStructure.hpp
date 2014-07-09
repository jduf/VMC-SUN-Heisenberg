#ifndef DEF_BANDSTRUCTURE
#define DEF_BANDSTRUCTURE

#include "Gnuplot.hpp"
#include "Lapack.hpp"

template<typename Type>
class BandStructure{
	public:
		/*!Constructor for a 1D BandStructure*/
		BandStructure(Matrix<Type> const& T, unsigned int const& L, unsigned int const& spuc, int const& bc);
		/*!Constructor for a 2D BandStructure*/
		BandStructure(Matrix<Type> const& T, unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, int const& bc);
		/*!Destructor*/
		virtual ~BandStructure(){}

	protected:
		Matrix<Type> Px_;//!< translation operator along x-axis
		Matrix<Type> Py_;//!< translation operator along y-axis
		Matrix<Type> TP_;//!< T+alpha*Px_+beta*Py_
		Matrix<std::complex<double> > evec_;//!< eigenvectors of TP_
		Vector<std::complex<double> > eval_;//!< eigenvalues of TP_
		unsigned int const Lx_;//!< number of unit cell along the x-axis
		unsigned int const Ly_;//!< number of unit cell along the y-axis
		unsigned int const spuc_;//!< site per unit cell
		int const bc_;//!< boundary condition

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int idx);
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		void compute_TP();

	private:
		/*!Forbids copy*/
		BandStructure(BandStructure const& bs);
		/*!Forbids assignment*/
		BandStructure& operator=(BandStructure bs);
};
	
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> const& T, unsigned int const& L, unsigned int const& spuc, int const& bc):
	Px_(T.row(),T.row(),0),
	Py_(T.row(),T.row(),0),
	TP_(T),
	Lx_(L),
	Ly_(0),
	spuc_(spuc),
	bc_(bc)
{
	compute_TP();
	//std::cout<<T*Px_-Px_*T<<std::endl;

	Lapack<Type> ES(TP_,false,'G');
	ES.eigensystem(eval_,&evec_);

	IOFiles w("spectrum.dat",true);
	for(unsigned int i(0);i<T.row();i++){
		w<<log(projection(Px_,i)).imag()<<" "<<projection(T,i).real()<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="plot 'spectrum.dat' u 1:2";
	gp.save_file();
}
	
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> const& T, unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, int const& bc):
	Px_(T.row(),T.row(),0),
	Py_(T.row(),T.row(),0),
	TP_(T),
	Lx_(Lx),
	Ly_(Ly),
	spuc_(spuc),
	bc_(bc)
{
	compute_TP();

	//std::cout<<T*Px_-Px_*T<<std::endl;
	//std::cout<<T*Py_-Py_*T<<std::endl;

	Lapack<Type> ES(TP_,false,'G');
	ES.eigensystem(eval_,&evec_);
	IOFiles w("spectrum.dat",true);
	for(unsigned int i(0);i<T.row();i++){
		w<<log(projection(Px_,i)).imag()<<" "<<log(projection(Py_,i)).imag()<<" "<<projection(T,i).real()<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="splot 'spectrum.dat' u 1:2:3";
	gp.save_file();
}

template<typename Type>
void BandStructure<Type>::compute_TP(){
	unsigned int tmp;
	if(Ly_){
		for(unsigned int j(0);j<Ly_;j++){
			for(unsigned int i(0);i<Lx_-1;i++){
				tmp = spuc_*(i + j*Lx_);
				for(unsigned int k(0);k<spuc_;k++){
					Px_(tmp+k, tmp+k+spuc_) = 1.0;
				}
			}
			tmp = spuc_*((Lx_-1) + j*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				Px_(tmp+k,spuc_*j*Lx_ + k) = bc_;
			}
		}
		Px_ *= Type(3.0);
		TP_ += Px_;
		for(unsigned int i(0);i<Lx_;i++){
			for(unsigned int j(0);j<Ly_-1;j++){
				tmp = spuc_*(i + j*Lx_);
				for(unsigned int k(0);k<spuc_;k++){
					Py_(tmp+k, tmp+spuc_*Lx_+k) = 1.0;
				}
			}
			tmp = spuc_*(i + (Ly_-1)*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				Py_(tmp+k, spuc_*i+k) = bc_;
			}
		}
		Py_ *= Type(7.0);
		TP_ += Py_;
	} else {
		for(unsigned int i(0); i<Lx_-1; i++){ 
			tmp = spuc_*i;
			for(unsigned int k(0);k<spuc_;k++){
				Px_(tmp+k,tmp+k+spuc_) = 1.0; 
			}
		}
		tmp = spuc_*(Lx_-1);
		for(unsigned int k(0);k<spuc_;k++){ Px_(tmp+k,k) = bc_; }
		TP_ += Px_;
	}
}

template<typename Type>
std::complex<double> BandStructure<Type>::projection(Matrix<Type> const& O, unsigned int idx){
	Vector<std::complex<double> > tmp(O.row(),0.0);
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		for(unsigned int j(0);j<O.col();j++){
			tmp(i) += O(i,j)*evec_(j,idx);
		}
	}
	for(unsigned int i(0);i<O.row();i++){
		out += std::conj(evec_(i,idx))*tmp(i);
	}
	return out;
}
#endif
