#ifndef DEF_SYSTEM2D
#define DEF_SYSTEM2D

#include "GenericSystem.hpp"

template<typename Type>
class System2D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System2D(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, unsigned int const& z, std::string const& filename, unsigned int const& sel0, unsigned int const& sel1);
		/*!Destructor*/
		virtual ~System2D()=0;

	protected:
		Matrix<Type> H_;			//!< matrix used to get the band structure
		unsigned int const Lx_;		//!< number of unit cell along the x-axis
		unsigned int const Ly_;		//!< number of unit cell along the y-axis
		Vector<unsigned int>* select_;
		Matrix<std::complex<double> > evec_;//!< eigenvectors of H+Tx+Ty

		void compute_TxTy();
		void compute_band_structure();
		void plot_band_structure();
		void select_eigenvectors();

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);

		Matrix<Type> Tx_;			//!< translation operator along x-axis
		Matrix<Type> Ty_;			//!< translation operator along y-axis
		Vector<double> e_;			//!< eigenvalue of the Hamiltonian T_
		Vector<double> px_;			//!< eigenvalue of the translation along x
		Vector<double> py_;			//!< eigenvalue of the translation along y
	private:
		unsigned int sel_[2];

};
	
/*{constructors*/
template<typename Type>
System2D<Type>::System2D(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, unsigned int const& z, std::string const& filename, unsigned int const& sel0, unsigned int const& sel1):
	GenericSystem<Type>(spuc,z,filename),
	Lx_(sqrt(Lx*this->n_/(Ly*spuc))),
	Ly_(sqrt(Ly*this->n_/(Lx*spuc))),
	select_(new Vector<unsigned int>[this->N_])
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
	Tx_.set(this->n_,this->n_,0);
	Ty_.set(this->n_,this->n_,0);
	unsigned int tmp;
	double t(1);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			tmp = this->spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<this->spuc_;k++){
				Tx_(tmp+k, tmp+k+this->spuc_) = t;
			}
		}
		tmp = this->spuc_*((Lx_-1) + j*Lx_);
		for(unsigned int k(0);k<this->spuc_;k++){
			Tx_(tmp+k,this->spuc_*j*Lx_ + k) = this->bc_*t;
		}
	}
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_-1;j++){
			tmp = this->spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<this->spuc_;k++){
				Ty_(tmp+k, tmp+this->spuc_*Lx_+k) = t;
			}
		}
		tmp = this->spuc_*(i + (Ly_-1)*Lx_);
		for(unsigned int k(0);k<this->spuc_;k++){
			Ty_(tmp+k, this->spuc_*i+k) = this->bc_*t;
		}
	}
}

template<typename Type>
void System2D<Type>::compute_band_structure(){
	e_.set(this->n_);
	Rand rnd(1e4);
	Matrix<Type> M;
	unsigned int iter(0);
	Matrix<Type> R(this->H_.row(),this->H_.col(),0);
	Matrix<int> nb;
	unsigned int s;
	unsigned int rs;
	int bc;
	for(unsigned int j(0); j<this->Ly_;j++){
		for(unsigned int i(0); i<this->Lx_;i++){
			bc = 1;
			s = this->spuc_*(i + j*this->Lx_);

			nb = this->get_neighbourg(s);
			rs = nb(0,0);
			bc*= nb(0,1);
			nb = this->get_neighbourg(rs);
			rs = nb(1,0);
			bc*= nb(1,1);
			for(unsigned int k(0); k<2*i;k++){
				nb = this->get_neighbourg(rs);
				rs = nb(1,0);
				bc*= nb(1,1);
			}
			if(j!=0){
				nb = this->get_neighbourg(rs);
				rs = nb(1,0);
				bc*= nb(1,1);
				for(unsigned int k(0); k<2*j-1;k++){
					nb = this->get_neighbourg(rs);
					rs = nb(2,0);
					bc*= nb(2,1);
				}
				nb = this->get_neighbourg(rs);
				rs = nb(3,0);
				bc*= nb(3,1);
			}
			nb = this->get_neighbourg(rs);
			R(s,rs) = bc;
			R(s+1,nb(0,0)) = bc*nb(0,1);
			R(s+2,nb(1,0)) = bc*nb(1,1);
			//std::cout<<s<<" "<<rs<<" "<<bc<<std::endl;
			//std::cout<<s+1<<" "<<nb(0,0)<<" "<<bc*nb(0,1)<<std::endl;
			//std::cout<<s+2<<" "<<nb(1,0)<<" "<<bc*nb(1,1)<<std::endl;
		}
	}

	do{
		/*may be optimized : */
		M  = H_;
		M += Tx_*Type(rnd.get()*rnd.get(1e6));
		M += Ty_*Type(rnd.get()*rnd.get(1e6));
		M += R*Type(rnd.get()*rnd.get(1e6));

		Vector<std::complex<double> > eval;
		Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);

		this->degenerate_ = false;
		for(unsigned int i(0);i<this->n_;i++){
			for(unsigned int j(i+1);j<this->n_;j++){
				if(are_equal(eval(i),eval(j),1e-10,1e-10)){
					std::cout<<eval(i)<<" "<<eval(j)<<std::endl;
					this->degenerate_ = true;
					i=j=this->n_;
				}
			}
		}
	} while ( this->degenerate_ && iter++<1e2);
	if( iter>=1e2 ){ std::cerr<<"void System2D<Type>::compute_band_structure() : degenerate"<<std::endl; }
	else {
		std::cout<<iter<<std::endl;
		Vector<unsigned int> index;
		Matrix<std::complex<double> > evec_tmp(evec_);
		for(unsigned int i(0);i<this->n_;i++){
			e_(i) = projection(H_,i).real();
		}
		e_.sort(std::less_equal<double>(),index);
		px_.set(this->n_);
		py_.set(this->n_);
		for(unsigned int i(0);i<this->n_;i++){
			for(unsigned int j(0);j<this->n_;j++){
				std::swap(evec_(i,j),evec_tmp(i,index(j)));
			}
		}
		for(unsigned int i(0);i<this->n_;i++){
			px_(i) = log(projection(Tx_,i)).imag();
			py_(i) = log(projection(Ty_,i)).imag();
		}
	}
}

template<typename Type>
void System2D<Type>::plot_band_structure(){
	IOFiles spectrum("spectrum.dat",true);
	for(unsigned int i(0);i<this->n_;i++){
		spectrum<<(are_equal(std::abs(px_(i)),M_PI,1e-12)?-M_PI:px_(i))<<" "<<(are_equal(std::abs(py_(i)),M_PI,1e-12)?-M_PI:py_(i))<<" "<<e_(i)<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="splot 'spectrum.dat' u 1:2:3";
	gp.save_file();
}

//template<typename Type>
//void System2D<Type>::select_eigenvectors(){
	////std::cout<<e_<<std::endl;
	//unsigned int iter(0);
	//double Px(0.0);
	//double Py(0.0);
	//for(unsigned int c(0);c<this->N_;c++){
		//unsigned int a(this->M_(c)-1);
		//unsigned int b(this->M_(c)-1);
		////std::cout<<a<<" "<<b<<std::endl;
		//do{b++;} while (b+1<this->n_ && are_equal(e_(b),e_(b-1)));
		////std::cout<<a<<" "<<b<<std::endl;
		//if(b!=this->M_(c)){ while(a>0 && are_equal(e_(a-1),e_(a))){a--;} }
		////std::cout<<a<<" "<<b<<std::endl;
		////std::cout<<e_.range(a,b)<<std::endl;
		//Vector<unsigned int> cnk;
		//Combination cbn;
		//cbn.set(this->M_(c)-a,b-a,cnk);
		//select_[c].set(this->M_(c));
		//for(unsigned int i(0);i<a;i++){ select_[c](i) = i; }
		////do{ 
		////for(unsigned int i(0);i+a<M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		////double Px(0.0);
		////double Py(0.0);
		////double  e(0.0);
		////for(unsigned int i(0);i<M_(c);i++){
		//////Px+= (are_equal(std::abs(px_(a+comb(j))),M_PI,1e-8)?0:px_(a+comb(j)));
		//////Py+= (are_equal(std::abs(py_(a+comb(j))),M_PI,1e-8)?0:py_(a+comb(j)));
		////e += e_(select_[c](i));
		////}
		////std::cout<<Px<<" "<<Py<<" "<<e<<std::endl;
		////} while (cbn.next());
		//while(iter++<sel_[c] && cbn.next());
		//for(unsigned int i(0);i+a<this->M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		//double  e(0.0);
		//for(unsigned int i(0);i<this->M_(c);i++){
			//Px+= are_equal(std::abs(px_(select_[c](i))),M_PI,1e-12)?0:px_(select_[c](i));
			//Py+= are_equal(std::abs(py_(select_[c](i))),M_PI,1e-12)?0:py_(select_[c](i));
			//e +=  e_(select_[c](i));
		//}
	//}
	//if(are_equal(Px,0.) && are_equal(Py,0.)){
		//std::cout<<iter<<std::endl;
	//}
//}

template<typename Type>
void System2D<Type>::select_eigenvectors(){
	//std::cout<<e_<<std::endl;
	unsigned int iter(0);
	for(unsigned int c(0);c<this->N_;c++){
		unsigned int a(this->M_(c)-1);
		unsigned int b(this->M_(c)-1);
		do{b++;} while (b+1<this->n_ && are_equal(e_(b),e_(b-1)));
		if(b!=this->M_(c)){ while(a>0 && are_equal(e_(a-1),e_(a))){a--;} }
		Vector<unsigned int> cnk;
		Combination cbn;
		cbn.set(this->M_(c)-a,b-a,cnk);
		select_[c].set(this->M_(c));
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
		for(unsigned int i(0);i+a<this->M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		double  e(0.0);
		double Px(0.0);
		double Py(0.0);
		for(unsigned int i(0);i<this->M_(c);i++){
			Px+= are_equal(std::abs(px_(select_[c](i))),M_PI,1e-12)?0:px_(select_[c](i));
			Py+= are_equal(std::abs(py_(select_[c](i))),M_PI,1e-12)?0:py_(select_[c](i));
			e +=  e_(select_[c](i));
		}
		if(!are_equal(Px,0.) || !are_equal(Py,0.)){ this->degenerate_ = true;}
		else {
			Px = 0.0;
			Py = 0.0;
			for(unsigned int i(a);i<this->M_(c);i++){
				Px+= are_equal(std::abs(px_(select_[c](i))),M_PI,1e-12)?0:px_(select_[c](i));
				Py+= are_equal(std::abs(py_(select_[c](i))),M_PI,1e-12)?0:py_(select_[c](i));
			}
			if(!are_equal(Px,0.) || !are_equal(Py,0.)){ this->degenerate_ = true;}
		}
	}
	//if(are_equal(Px,0.) && are_equal(Py,0.)){
	//this->degenerate_ = false;
	//for(unsigned int c(0);c<this->N_;c++){
	//std::cout<<px_(select_[c](this->M_(c)-1))<<" "<<py_(select_[c](this->M_(c)-1))<<" ";
	//}
	//} else {
	//std::cout<<sel_[0]<<" "<<sel_[1]<<std::endl;
	//this->degenerate_ = true;
	//}
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
