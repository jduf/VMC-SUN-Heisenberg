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
		Matrix<Type> H_;		//!< matrix used to get the band structure
		unsigned int const Lx_;	//!< number of unit cell along the x-axis
		unsigned int const Ly_;	//!< number of unit cell along the y-axis
		Matrix<std::complex<double> > evec_;//!< eigenvectors of H+Tx+Ty
		Vector<unsigned int>* select_;

		/*!Plot the band structure E(px,py)*/
		void plot_band_structure();
		/*!Create the selection of optimal eigenvectors*/
		void select_eigenvectors();

		void diagonalize(bool simple);

	private:
		Matrix<Type> Tx_;	//!< translation operator along x-axis
		Matrix<Type> Ty_;	//!< translation operator along y-axis
		Vector<double> px_;	//!< eigenvalue of Tx
		Vector<double> py_;	//!< eigenvalue of Ty
		Vector<double> e_;	//!< eigenvalue of H_

		/*!Compute the translation operators*/
		void compute_TxTy();
		/*!Diagonalize H_*/
		bool simple_diagonalization();
		/*!Diagonalize H_+T_ => compute the band structure E(p)*/
		bool full_diagonalization();
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);

		/*might be removed... used for magnetization... but not good solution*/
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
	if(this->n_==this->spuc_*Lx_*Ly_){
		this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		this->status_--;
	} else {
		std::cerr<<"System2D<Type> : the cluster is impossible, n must be a"<<std::endl; 
		std::cerr<<"               : multiple of "<<Lx*Ly*spuc<<" ("<<Lx<<"x"<<Ly<<"x"<<spuc<<")"<<std::endl; 
	}
}

template<typename Type>
System2D<Type>::~System2D(){
	if(select_){ delete[] select_;}
}
/*}*/

/*{protected methods*/
template<typename Type>
void System2D<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_--; } }
	else { if(full_diagonalization()){ this->status_--; } }
}

template<typename Type>
void System2D<Type>::plot_band_structure(){
	full_diagonalization();

	IOFiles spectrum("spectrum.dat",true);
	for(unsigned int i(0);i<this->n_;i++){
		spectrum<<(are_equal(std::abs(px_(i)),M_PI,1e-12)?-M_PI:px_(i))<<" "<<(are_equal(std::abs(py_(i)),M_PI,1e-12)?-M_PI:py_(i))<<" "<<e_(i)<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.range("x","-pi","pi");
	gp.range("y","-pi","pi");
	gp.range("z","-5","5");
	gp+="splot 'spectrum.dat' u 1:2:3";
	gp.save_file();
}

template<typename Type>
void System2D<Type>::select_eigenvectors(){
/*{*/
	/*!
	  unsigned int iter(0);
	  Matrix<double> pxpy[2];
	  this->degenerate_ = false; //normally useless
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
	  while(iter++<sel_[c] && cbn.next());
	  for(unsigned int i(0);i+a<this->M_(c);i++){ select_[c](a+i) = a+cnk(i); }
	  double  e(0.0);
	  double Px(0.0);
	  double Py(0.0);
	  for(unsigned int i(0);i<this->M_(c);i++){
	  Px+= are_equal(std::abs(px_(select_[c](i))),M_PI,1e-12,1e-12)?0:px_(select_[c](i));
	  Py+= are_equal(std::abs(py_(select_[c](i))),M_PI,1e-12,1e-12)?0:py_(select_[c](i));
	  e += e_(select_[c](i));
	  }
	  if(!are_equal(Px,0.,1e-14) || !are_equal(Py,0.,1e-14)){ 
	  this->degenerate_ = true;
	  } else {
	  pxpy[c].set(this->M_(c)-a,2);
	  for(unsigned int i(a);i<this->M_(c);i++){
	  pxpy[c](i-a,0) = px_(select_[c](i));
	  pxpy[c](i-a,1) = py_(select_[c](i));
	  }
	  }
	  }
	  if(!this->degenerate_){
	//for(unsigned int c(0);c<this->N_;c++){
	//for(unsigned int i(0);i<pxpy[c].row();i++){
	//std::cout<<pxpy[c](i,0)<<" "<<pxpy[c](i,1)<<" ";
	//}
	//}
	}

*/
	/*}*/
	unsigned int c(0);
	unsigned int a(this->M_(c)-1);
	unsigned int b(this->M_(c)-1);
	do{b++;} while (b+1<this->n_ && are_equal(e_(b),e_(b-1)));
	if(b!=this->M_(c)){ while(a>0 && are_equal(e_(a-1),e_(a))){a--;} }
	std::cout<<a<<" "<<b<<std::endl;
	std::cout<<e_<<std::endl;
	for(unsigned int i(a-1);i<b+1;i++){
		std::cout<<i<<" "<<chop(e_(i))<<" "<<px_(i)<<" "<<py_(i)<<std::endl;
	}
	Matrix<unsigned int> pair(b-a,2,0);
	for(unsigned int i(a);i<b;i++){
		for(unsigned int j(a);j<b;j++){
			if(are_equal(e_(i),e_(j)) && are_equal(px_(i),-px_(j)) && are_equal(py_(i),-py_(j))){
				pair(i-a,0) = i;
				pair(i-a,1) = j;
			}
		}
	}
	std::cout<<pair<<std::endl;
	for(unsigned int j(0);j<pair.row()/2;j++){
		double n1(0);
		double n2(0);
		std::complex<double> tmp1;
		std::complex<double> tmp2;
		for(unsigned int i(0);i<this->n_;i++){
			tmp1 = evec_(i,pair(j,0)) + evec_(i,pair(j,1));
			tmp2 = evec_(i,pair(j,0)) - evec_(i,pair(j,1));
			evec_(i,pair(j,0)) = tmp1;
			evec_(i,pair(j,1)) = tmp2;
			n1 += norm_squared(tmp1);
			n2 += norm_squared(tmp2);
		}
		for(unsigned int i(0);i<this->n_;i++){
			evec_(i,pair(j,0)) /= sqrt(n1);
			evec_(i,pair(j,1)) /= sqrt(n2);
		}
	}
	//std::cout<<(evec_*evec_.trans_conj()).chop().diag()<<std::endl;
	std::cout<<"this selection method is not enough general !"<<std::endl;
}
/*}*/

/*{private methods*/
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
bool System2D<Type>::simple_diagonalization(){
	Vector<double> eval;
	Lapack<Type>(H_,false,(this->ref_(1)==1?'S':'H')).eigensystem(eval,true);
	for(unsigned int c(0);c<this->N_;c++){
		if(are_equal(eval(this->M_(c)),eval(this->M_(c)-1),1e-12)){
			std::cerr<<"bool System2D<Type>::simple_diagonalization() :"
				" degenerate at the Fermi level"<<std::endl;
			return false;
		}
	}
	return true;
}

template<typename Type>
bool System2D<Type>::full_diagonalization(){
	compute_TxTy();
	Matrix<Type> M(H_);
	M += Tx_*Type(3.0);
	M += Ty_*Type(7.0);
	Vector<std::complex<double> > eval;
	Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(are_equal(eval(i),eval(j),1e-10,1e-10)){
				std::cerr<<"bool System2D<Type>::full_diagonalization() :"
					"eigenvalue "<<i<<" and "<<j<<" degenerate"<<std::endl;
				return false;
			}
		}
	}
	Vector<unsigned int> index;
	e_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){ e_(i) = projection(H_,i).real(); }
	e_.sort(std::less_equal<double>(),index);

	Matrix<std::complex<double> > evec_tmp(evec_);
	Vector<std::complex<double> > eval_tmp(eval);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			std::swap(evec_(i,j),evec_tmp(i,index(j)));
		}
		std::swap(eval(i),eval_tmp(index(i)));
	}

	px_.set(this->n_);
	py_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		px_(i) = log(projection(Tx_,i)).imag();
		py_(i) = log(projection(Ty_,i)).imag();
	}
	return true;
}

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
