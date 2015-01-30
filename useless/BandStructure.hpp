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

		void diagonalize_subspace_Tx(unsigned int const& a, unsigned int const&b, Matrix<std::complex<double> >& evec, Vector<double> const& eval);
		void diagonalize_subspace_Ty(unsigned int const& a, unsigned int const&b, Matrix<std::complex<double> >& evec, Vector<double> const& eval);
		void diagonalize_everything(Matrix<std::complex<double> >& evec, Matrix<std::complex<double> >& eval);
		void check(Matrix<std::complex<double> >O, Matrix<std::complex<double> > const& evec, Matrix<std::complex<double> > const& eval, unsigned int idx);
		void save();

		void compute_band_structure();
		Matrix<std::complex<double> > get_evec() const { return evec_; };
		Vector<double> get_E() const { return E_; };
		Vector<double> get_px() const { return px_; };
		Vector<double> get_py() const { return py_; };

	protected:
		Matrix<Type> mat_;			//!< matrix used to get the band structure
		Matrix<Type> Tx_;			//!< translation operator along x-axis
		Matrix<Type> Ty_;			//!< translation operator along y-axis
		Matrix<Type> evec_;			//!< translation operator along y-axis
		Vector<double> E_;			//!< eigenvalue of the Hamiltonian T_
		Vector<double> px_;			//!< eigenvalue of the translation along x
		Vector<double> py_;			//!< eigenvalue of the translation along y
		unsigned int const Lx_;		//!< number of unit cell along the x-axis
		unsigned int const Ly_;		//!< number of unit cell along the y-axis
		unsigned int const spuc_;	//!< site per unit cell
		int const bc_;				//!< boundary condition
		unsigned int const dim_;	//!< dimension of the lattice

	private:
		/*!Forbids copy*/
		BandStructure(BandStructure const& bs);
		/*!Forbids assignment*/
		BandStructure& operator=(BandStructure bs);

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, Matrix<Type> const& evec, unsigned int const& idx, unsigned int const& jdx) const;
		/*!Creates Tx and Ty*/
		void init();
};
	
/*{constructors*/
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
{ init(); 

}
	
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
{ init(); }

template<typename Type>
void BandStructure<Type>::init(){
	unsigned int tmp;
	double t;
	switch(dim_){
		case 1:
			{
				t = 3.0;
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
/*}*/

/*{public methods*/
template<typename Type>
void BandStructure<Type>::compute_band_structure(){
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > M((dim_==1)?mat_+Tx_:mat_+Tx_+Ty_);
	Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);
	for(unsigned int i(0);i<mat_.row();i++){
		E_(i) = projection(mat_,evec_,i,i).real();
		px_(i) = log(projection(Tx_,evec_,i,i)).imag();
		if(dim_==2){ py_(i) = log(projection(Ty_,evec_,i,i)).imag(); }
	}
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
void BandStructure<Type>::diagonalize_subspace_Tx(unsigned int const& a, unsigned int const&b, Matrix<std::complex<double> >& evec, Vector<double> const& eval){
	if(b-a!=1){
		Matrix<std::complex<double> > evec_old(evec.row(),b-a);
		for(unsigned int i(0);i<evec.row();i++){
			for(unsigned int j(0);j<b-a;j++){
				evec_old(i,j) = evec(i,a+j);
			}
		}

		Matrix<std::complex<double> > tmp_mat(b-a,b-a);
		for(unsigned int i(0);i<b-a;i++){
			for(unsigned int j(0);j<b-a;j++){
				tmp_mat(i,j) = projection(Tx_,evec_old,i,j);
			}
		}

		Lapack<std::complex<double> > tmp_diag(tmp_mat,false,'G');
		Vector<std::complex<double> > tmp_eval;
		Matrix<std::complex<double> > tmp_evec;
		tmp_diag.eigensystem(tmp_eval,&tmp_evec);

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

		for(unsigned int i(0);i<b-a;i++){
			px_(a+i)= log(tmp_eval(i)).imag();
			E_(a+i) = eval(a+i);
		}
	} else {
		px_(a)= log(projection(Tx_,evec,a,a)).imag();
		E_(a) = eval(a);
	}
}

//template<typename Type>
//void BandStructure<Type>::diagonalize_everything(Matrix<std::complex<double> >& evec, Matrix<std::complex<double> >& eval){
	//Matrix<Type>* O(new Matrix<Type>[dim_]);
	//O[0] = Tx_;
	//if(dim_==2){ O[1] = Ty_; }
//
	//unsigned int a;
	//unsigned int b;
	//unsigned int n(evec.row());//can be replaced by n_
	//for(unsigned int o(0);o<dim_;o++){
//
			//std::cout<<"mat "<<(evec.transpose()*evec).chop()<<std::endl<<std::endl;
		//a = 0;
		//b = 0;
		//while(a != n){
			//do{b++;}
			//while(b != n && are_equal(eval(b,o),eval(b-1,o),1e-14));
			//std::cout<<o<<" "<<a<<" "<<b<<std::endl;
//
			//if(b-a>1){
				///*copy the old eigenvectors*/
				//Matrix<std::complex<double> > evec_old(evec.row(),b-a);
				//for(unsigned int i(0);i<evec.row();i++){
					//for(unsigned int j(0);j<b-a;j++){
						//evec_old(i,j) = evec(i,a+j);
					//}
				//}
				///*project the old eigenvectors on the operator O[o]*/
				//Matrix<std::complex<double> > tmp_mat(b-a,b-a);
				//for(unsigned int i(0);i<b-a;i++){
					//for(unsigned int j(0);j<b-a;j++){
						//tmp_mat(i,j) = projection(O[o],evec_old,i,j);
					//}
				//}
				///*diagonalize the subspace and compute the eigenvectors*/
				//Lapack<std::complex<double> > tmp_diag(tmp_mat,false,'G');
				//Matrix<std::complex<double> > tmp_evec;
				//Vector<std::complex<double> > tmp_eval;
				//tmp_diag.eigensystem(tmp_eval,&tmp_evec);
				//std::cout<<tmp_eval<<std::endl;
				///*copy the eigenvalue in a sorted order*/
				//for(unsigned int i(0);i<b-a;i++){ 
					//tmp_eval(i) = log(tmp_eval(i));
				//}
				//Vector<unsigned int> index;
				//tmp_eval.sort(&less_equal_complex, index);
				//for(unsigned int i(0);i<b-a;i++){ 
					//eval(a+i,o+1) = tmp_eval(i); 
				//}
				///*create the linear combination of the old vectors, the new
				// * vectors are eigenvectors of the operator O[o]*/
				//for(unsigned int i(0);i<evec.row();i++){
					//for(unsigned int j(0);j<b-a;j++){
						//evec(i,a+j) = 0.0;
						//for(unsigned int k(0);k<b-a;k++){
							//evec(i,a+j)	+= evec_old(i,k)*tmp_evec(k,index(j));
						//}
					//}
				//}
				//Matrix<std::complex<double> > evecinv(evec);
				//Lapack<std::complex<double> > (evecinv,false,'G').inv();
				//std::cout<<(evecinv*O[o]*evec).chop()<<std::endl<<std::endl;
//
			//} else {
				//eval(a,o+1) = log(projection(O[o],evec,a,a));
			//}
			//a=b;
		//}
	//}
	//delete[] O;
//}


inline bool sort_complex(std::complex<double> a, std::complex<double> b){
	if(are_equal(a.real(),b.real(),1e-12)){ return a.imag()<=b.imag(); }
	else{ return a.real()<=b.real(); }
}

template<typename Type>
void BandStructure<Type>::diagonalize_everything(Matrix<std::complex<double> >& evec, Matrix<std::complex<double> >& eval){
	Matrix<Type>* O(new Matrix<Type>[dim_]);
	O[0] = Tx_;
	if(dim_==2){ O[1] = Ty_; }

	unsigned int a;
	unsigned int b;
	unsigned int n(evec.row());//can be replaced by n_
	for(unsigned int o(0);o<dim_;o++){
		a = 0;
		b = 0;
		while(a != n){
			do{b++;}
			while(b != n && are_equal(eval(b,o),eval(b-1,o),1e-14));
			std::cout<<o<<" "<<b-a<<std::endl;

			if(b-a>1){
				/*copy the old eigenvectors*/
				Matrix<std::complex<double> > evec_old(evec.row(),b-a);
				for(unsigned int i(0);i<evec.row();i++){
					for(unsigned int j(0);j<b-a;j++){
						evec_old(i,j) = evec(i,a+j);
					}
				}
				/*project the old eigenvectors on the operator O[o]*/
				Matrix<std::complex<double> > tmp_mat(b-a,b-a);
				for(unsigned int i(0);i<b-a;i++){
					for(unsigned int j(0);j<b-a;j++){
						tmp_mat(i,j) = projection(O[o],evec_old,i,j);
					}
				}
				/*diagonalize the subspace and compute the eigenvectors*/
				Vector<std::complex<double> > tmp_eval;
				Matrix<std::complex<double> > r_tmp_evec;
				Matrix<std::complex<double> > l_tmp_evec;
				Lapack<std::complex<double> >(tmp_mat,false,'G').eigensystem(tmp_eval,&r_tmp_evec,&l_tmp_evec);
				//std::cout<<tmp_eval<<std::endl;
				/*copy the eigenvalue in a sorted order*/
				//for(unsigned int i(0);i<b-a;i++){ 
					//tmp_eval(i) = log(tmp_eval(i));
				//}
				Vector<unsigned int> index;
				tmp_eval.sort(&sort_complex,index);
				for(unsigned int i(0);i<b-a;i++){ 
					eval(a+i,o+1) = tmp_eval(index(i)); 
				}
				/*create the linear combination of the old vectors, the new
				 * vectors are eigenvectors of the operator O[o]*/
				for(unsigned int i(0);i<evec.row();i++){
					for(unsigned int j(0);j<b-a;j++){
						evec(i,a+j) = 0.0;
						for(unsigned int k(0);k<b-a;k++){
							evec(i,a+j)	+= evec_old(i,k)*r_tmp_evec(k,index(j));
						}
					}
				}
				Matrix<std::complex<double> > evecinv(evec);
				Lapack<std::complex<double> >(evecinv,false,'G').inv();
				std::cout<<(evecinv*O[o]*evec).chop()<<std::endl<<std::endl;
			} else {
				eval(a,o+1) = log(projection(O[o],evec,a,a));
			}
			a=b;
		}
		//std::cout<<"here"<<std::endl;
		//for(unsigned int i(0);i<evec.row();i++){
			//for(unsigned int j(0);j<evec.row();j++){
				//std::complex<double> tmp(0);
				//for(unsigned int k(0);k<evec.row();k++){
					//tmp += O[o](i,k)*evec(k,j);
				//}
				//std::cout<<tmp-evec(i,j)*eval(i,o)<<std::endl;
			//}
		//}
	}
	std::cout<<eval<<std::endl;
	check(O[0],evec,eval,1);
	delete[] O;
}

template<typename Type>
void BandStructure<Type>::check(Matrix<std::complex<double> > O, Matrix<std::complex<double> > const& evec, Matrix<std::complex<double> > const& eval, unsigned int idx){
	for(unsigned int i(0);i<O.row();i++){
		O(i,i) -= eval(0,idx);
	}
	std::complex<double> tmp;
	for(unsigned int i(0);i<O.row();i++){
		for(unsigned int j(0);j<O.col();j++){
			tmp = 0;
			for(unsigned int k(0);k<O.col();k++){
				tmp += O(i,k)*evec(k,j);
			}
			//if(!are_equal(tmp,0)){ std::cout<<i<<" "<<j<<std::endl;}
		}
	}
	std::cout<<(O*evec).chop()<<std::endl;
}
/*}*/

/*{private methods*/
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
/*}*/
#endif
