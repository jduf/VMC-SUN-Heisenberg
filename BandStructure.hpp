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
		void dostuff(unsigned int const& a, unsigned int const&b, Matrix<std::complex<double> >& evec, Vector<double> const& eval);
		Vector<double> check(Matrix<std::complex<double> >& evec);
		void save();

	protected:
		Matrix<Type> mat_;			//!< matrix used to get the band structure
		Matrix<Type> Tx_;			//!< translation operator along x-axis
		Matrix<Type> Ty_;			//!< translation operator along y-axis
		Vector<double> E_;			//!< eigenvalue of the Hamiltonian T_
		Vector<double> px_;			//!< eigenvalue of the translation along x
		Vector<double> py_;			//!< eigenvalue of the translation along y
		unsigned int const Lx_;		//!< number of unit cell along the x-axis
		unsigned int const Ly_;		//!< number of unit cell along the y-axis
		unsigned int const spuc_;	//!< site per unit cell
		int const bc_;				//!< boundary condition
		unsigned int const dim_;

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, Matrix<Type> const& evec, unsigned int const& idx, unsigned int const& jdx) const;
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		void compute_translation_operator();
		void compute_band_structure();

	private:
		/*!Forbids copy*/
		BandStructure(BandStructure const& bs);
		/*!Forbids assignment*/
		BandStructure& operator=(BandStructure bs);
};
	
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> const& T, unsigned int const& L, unsigned int const& spuc, int const& bc):
	mat_(T),
	Tx_(T.row(),T.row(),0),
	E_(T.row()),
	px_(T.row()),
	Lx_(L),
	Ly_(0),
	spuc_(spuc),
	bc_(bc),
	dim_(1)
{ compute_translation_operator(); }
	
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> const& T, unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, int const& bc):
	mat_(T),
	Tx_(T.row(),T.row(),0),
	Ty_(T.row(),T.row(),0),
	E_(T.row()),
	px_(T.row()),
	py_(T.row()),
	Lx_(Lx),
	Ly_(Ly),
	spuc_(spuc),
	bc_(bc),
	dim_(2)
{}

template<typename Type>
void BandStructure<Type>::compute_translation_operator(){
	unsigned int tmp;
	double t;
	switch(dim_){
		case 1:
			{
				t = 1.0;
				for(unsigned int i(0); i<Lx_-1; i++){ 
					tmp = spuc_*i;
					for(unsigned int k(0);k<spuc_;k++){
						Tx_(tmp+k,tmp+k+spuc_) = t; 
					}
				}
				tmp = spuc_*(Lx_-1);
				for(unsigned int k(0);k<spuc_;k++){ Tx_(tmp+k,k) = bc_*t; }
			}break;
		case 2:
			{
				t = 3.0;
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
				t = 7.0;
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
			}break;
		default:{ std::cerr<<"BandStructure : dim_ undefined"<<std::endl;}
	}
}

template<typename Type>
void BandStructure<Type>::compute_band_structure(){
	Lapack<Type> ES(((dim_==1)?mat_+Tx_:mat_+Tx_+Ty_),true,'G');
	Matrix<std::complex<double> > evec;
	Vector<std::complex<double> > eval;
	ES.eigensystem(eval,&evec);
	for(unsigned int i(0);i<mat_.row();i++){
		E_(i) = projection(mat_,evec,i,i).real();
		px_(i) = log(projection(Tx_,evec,i,i)).imag();
		if(dim_==2){ py_(i) = log(projection(Ty_,evec,i,i)).imag(); }
		//px_(i) = log(projection(Px_,i)).imag();
		//if(dim_==2){ py_(i) = log(projection(Ty_,i)).imag(); }
	}
}

template<typename Type>
std::complex<double> BandStructure<Type>::projection(Matrix<Type> const& O, Matrix<Type> const& evec, unsigned int const& idx, unsigned int const& jdx) const {
	std::complex<double> tmp;
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		tmp = 0.0;
		for(unsigned int j(0);j<O.col();j++){
			tmp += O(i,j)*evec(j,idx);
		}
		out += std::conj(evec(i,jdx))*tmp;
	}
	return out;
}

template<typename Type>
void BandStructure<Type>::save(){
	IOFiles spectrum("spectrum.dat",true);
	for(unsigned int i(0);i<mat_.row();i++){
		spectrum<<px_(i)<<" ";
		if(dim_==2){spectrum<<py_(i)<<" ";}
		spectrum<<E_(i)<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	if(dim_==2){ gp+="splot 'spectrum.dat' u 1:2:3";}
	else{ gp+="plot 'spectrum.dat' u 1:2";}
	gp.save_file();
}

template<typename Type>
void BandStructure<Type>::dostuff(unsigned int const& a, unsigned int const&b, Matrix<std::complex<double> >& evec, Vector<double> const& eval){
	if(b-a!=1){
		Matrix<std::complex<double> > tmp_mat(b-a,b-a);
		for(unsigned int i(0);i<b-a;i++){
			for(unsigned int j(0);j<b-a;j++){
				tmp_mat(i,j) = projection(Tx_,evec,a+i,a+j);
			}
		}
		Lapack<std::complex<double> > tmp_diag(tmp_mat,false,'G');
		Vector<std::complex<double> > tmp_eval;
		Matrix<std::complex<double> > tmp_evec;
		tmp_diag.eigensystem(tmp_eval,&tmp_evec);

		for(unsigned int i(0);i<b-a;i++){
			px_(a+i) = log(tmp_eval(i)).imag();
			E_(a+i) = eval(a+i);
		}
		
		Matrix<std::complex<double> > evec_old(evec.row(),b-a);

		for(unsigned int i(0);i<evec.row();i++){
			for(unsigned int j(0);j<b-a;j++){
				evec_old(i,j) = evec(i,a+j);
			}
		}
		std::complex<double> tmp;
		for(unsigned int i(0);i<evec.row();i++){
			for(unsigned int j(0);j<b-a;j++){
				tmp = 0.0;
				for(unsigned int k(0);k<b-a;k++){
					tmp += evec_old(i,k)*tmp_evec(k,j);
				}
				evec(i,a+j) = tmp;
			}
		}
	} else {
		px_(a) = log(projection(Tx_,evec,a,a)).imag();
		E_(a) = eval(a);
	}
}

template<typename Type>
Vector<double> BandStructure<Type>::check(Matrix<std::complex<double> >& evec){
	Matrix<std::complex<double> > m(evec.trans_conj()*Tx_*evec);
	//std::cout<<m.chop()<<std::endl;
	//std::cout<<(evec.trans_conj()*mat_*evec).diag().chop()<<std::endl;
	Vector<double> p(evec.row());
	for(unsigned int i(0);i<evec.row();i++){ p(i) = log(m(i,i)).imag(); }
	return p;
}
#endif
