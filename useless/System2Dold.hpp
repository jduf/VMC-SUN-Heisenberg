#ifndef DEF_SYSTEM2D
#define DEF_SYSTEM2D

#include "GenericSystem.hpp"

/*{Description*/
/*!The rules that the arguments of the constructor must obey are the folloing :
 *
 * + **the unit is the lattice spacing**
 * + **vector expressed in the Cartesian coordinates**
 *
 * These two rules must apply to :
 *
 * + LxLy
 * + ab
 * + linear_jump
 *
 * This convention allows the use of local basis in the definitions of :
 *
 * + virtual Vector<double> get_pos_in_lattice(unsigned int const& i) const = 0;
 * + virtual unsigned int unit_cell_index(Vector<double> const& x) const = 0;
 */
/*}*/
template<typename Type>
class System2D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System2D(Matrix<double> const& LxLy, Matrix<double> const& ab, Vector<double> const& linear_jump, unsigned int const& spuc, unsigned int const& z, std::string const& filename);
		/*!Destructor*/
		virtual ~System2D() = default;

	protected:
		unsigned int xloop_;		//!< linear size of the lattice (number jumps before coming back to the origin)
		Vector<double> linear_jump_;//!< set_pos_LxLy(xloop_*linear_jump_) must come back to (0,0)
		Matrix<Type> H_;			//!< matrix used to get the band structure
		Matrix<double> ab_;  		//!< basis of the unit cel
		Matrix<double> LxLy_;		//!< basis of the whole lattice
		Matrix<double> dir_nn_LxLy_;//!<direction of the nearest neighbour in the LxLy basis
		Matrix<std::complex<double> > evec_;//!< eigenvectors of H+Tx+Ty
		double const eq_prec_;		//!< precision for equality (important for matchinf position in lattice)

		/*!Plots the band structure E(px,py)*/
		void plot_band_structure();
		/*!Creates the selection of optimal eigenvectors*/
		void select_eigenvectors();

		void diagonalize(bool simple);

		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Returns the position of the site i in the lattice basis*/
		virtual Vector<double> get_pos_in_lattice(unsigned int const& i) const = 0;
		/*!Returns the index of the position x the unit cell basis (a,b)*/
		virtual unsigned int unit_cell_index(Vector<double> const& x) const = 0;
		/*!Returns the index of the site i in the unit cell*/
		unsigned int get_site_in_unit_cell(unsigned int const& i) const;
		/*!Makes a change of basis to express the position x in the lattice basis (Lx,Ly)*/
		void set_pos_LxLy(Vector<double>& x) const;
		/*!Makes a change of basis to express the position x in the unit cell (a,b)*/
		void set_pos_ab(Vector<double>& x) const;

	private:
		Matrix<Type> Tx_;		//!< translation operator along x-axis
		Matrix<Type> Ty_;		//!< translation operator along y-axis
		Vector<double> px_;		//!< eigenvalue of Tx
		Vector<double> py_;		//!< eigenvalue of Ty
		Vector<double> e_;		//!< eigenvalue of H_
		Matrix<double> inv_LxLy_;//!<inverse of the matrix LxLy_
		Matrix<double> inv_ab_;	//!< inverse of the matrix ab_

		/*!Computes the translation operators*/
		void compute_TxTy();
		/*!Diagonalizes H_*/
		bool simple_diagonalization();
		/*!Diagonalizes H_+T_ => compute the band structure E(p)*/
		bool full_diagonalization();
		/*!Evaluates the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);

		virtual Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const = 0;
		virtual void try_neighbourg(Vector<double>& tn, unsigned int const& i) const = 0;

		/*!Makes sure that 0<=x_i<1, i=1,2*/
		void set_pos(Vector<double>& x) const;
		/*!Reset x so that it belongs to the lattice (Lx,Ly)*/
		bool set_in_LxLy(Vector<double>& x) const;
};

/*{constructors*/
template<typename Type>
System2D<Type>::System2D(Matrix<double> const& LxLy, Matrix<double> const& ab, Vector<double> const& linear_jump, unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	xloop_(0),
	linear_jump_(linear_jump),
	ab_(ab),
	LxLy_(LxLy),
	eq_prec_(1e-12),
	inv_LxLy_(2,2),
	inv_ab_(2,2)
{
	if(LxLy_.size()){
		inv_LxLy_(0,0) = LxLy_(1,1);
		inv_LxLy_(1,0) =-LxLy_(1,0);
		inv_LxLy_(0,1) =-LxLy_(0,1);
		inv_LxLy_(1,1) = LxLy_(0,0);
		inv_LxLy_/=(LxLy_(0,0)*LxLy_(1,1)-LxLy_(1,0)*LxLy_(0,1));

		inv_ab_(0,0) = ab_(1,1);
		inv_ab_(1,0) =-ab_(1,0);
		inv_ab_(0,1) =-ab_(0,1);
		inv_ab_(1,1) = ab_(0,0);
		inv_ab_/=(ab_(0,0)*ab_(1,1)-ab_(1,0)*ab_(0,1));

		Vector<double> x(2);
		do{//might be a more efficient way
			xloop_++;
			x = linear_jump_*xloop_;
			set_pos_LxLy(x);
		} while(!my::are_equal(x(0),0) || !my::are_equal(x(1),0));

		Vector<double> Lx(2);
		Lx(0) = LxLy_(0,0);
		Lx(1) = LxLy_(1,0);
		set_pos_ab(Lx);
		Vector<double> Ly(2);
		Ly(0) = LxLy_(0,1);
		Ly(1) = LxLy_(1,1);
		set_pos_ab(Ly);

		Vector<double> zero(2,0);
		if( my::are_equal(Lx,zero) &&  my::are_equal(Ly,zero) ){
			if(this->spuc_){ this->status_--; }
			else { std::cerr<<__PRETTY_FUNCTION__<<" : the unit cell contains 0 site"<<std::endl; }
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : the unit cell doesn't fit into the cluster"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the number of site doesn't fit into the cluster"<<std::endl; }
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
		spectrum<<(my::are_equal(std::abs(px_(i)),M_PI,1e-12)?-M_PI:px_(i))<<" "<<(my::are_equal(std::abs(py_(i)),M_PI,1e-12)?-M_PI:py_(i))<<" "<<e_(i)<<IOFiles::endl;
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
	  do{ b++; } while (b+1<this->n_ && my::are_equal(e_(b),e_(b-1)));
	  if(b!=this->M_(c)){ while(a>0 && my::are_equal(e_(a-1),e_(a))){ a--; } }
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
	  Px+= my::are_equal(std::abs(px_(select_[c](i))),M_PI,1e-12,1e-12)?0:px_(select_[c](i));
	  Py+= my::are_equal(std::abs(py_(select_[c](i))),M_PI,1e-12,1e-12)?0:py_(select_[c](i));
	  e += e_(select_[c](i));
	  }
	  if(!my::are_equal(Px,0.,1e-14) || !my::are_equal(Py,0.,1e-14)){
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
	do{ b++; } while (b+1<this->n_ && my::are_equal(e_(b),e_(b-1)));
	if(b!=this->M_(c)){ while(a>0 && my::are_equal(e_(a-1),e_(a))){ a--; } }
	std::cout<<a<<" "<<b<<std::endl;
	std::cout<<e_<<std::endl;
	for(unsigned int i(a-1);i<b+1;i++){
		std::cout<<i<<" "<<my::chop(e_(i))<<" "<<px_(i)<<" "<<py_(i)<<std::endl;
	}
	Matrix<unsigned int> pair(b-a,2,0);
	for(unsigned int i(a);i<b;i++){
		for(unsigned int j(a);j<b;j++){
			if(my::are_equal(e_(i),e_(j)) && my::are_equal(px_(i),-px_(j)) && my::are_equal(py_(i),-py_(j))){
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
			n1 += my::norm_squared(tmp1);
			n2 += my::norm_squared(tmp2);
		}
		for(unsigned int i(0);i<this->n_;i++){
			evec_(i,pair(j,0)) /= sqrt(n1);
			evec_(i,pair(j,1)) /= sqrt(n2);
		}
	}
	//std::cout<<(evec_*evec_.trans_conj()).chop().diag()<<std::endl;
	std::cout<<"this selection method is not enough general !"<<std::endl;
}

template<typename Type>
Matrix<int> System2D<Type>::get_neighbourg(unsigned int const& i) const {
	/*!tn is the position in the lattice of site i*/
	Vector<double> tn(get_pos_in_lattice(i));
	set_pos_LxLy(tn);
	set_in_LxLy(tn);

	/*!nn* are the nearest neighbours */
	Vector<double>* nn(new Vector<double>[this->z_]);
	Matrix<int> nb(this->z_,2);
	std::vector<unsigned int> dir(this->z_);
	for(unsigned int d(0);d<this->z_;d++){
		dir[d]= d;
		nn[d] = tn;
		nn[d]+= vector_towards(i,d);
		nb(d,1) = set_in_LxLy(nn[d]);
	}

	/*!tn will be the trial neighbour*/
	tn.set(2,0);
	unsigned int j(0);
	while(dir.size() && j<this->n_+2){
		for(unsigned int d(0);d<dir.size();d++){
			if(my::are_equal(tn,nn[dir[d]])){
				nb(dir[d],0) = j;
				dir.erase(dir.begin()+d);
			}
		}
		j++;
		/*!go to the trial next neighbour*/
		try_neighbourg(tn,j);
		set_in_LxLy(tn);
	}
	assert(j<this->n_+1);
	delete[] nn;
	return nb;
}

template<typename Type>
unsigned int System2D<Type>::get_site_in_unit_cell(unsigned int const& i) const {
	Vector<double> x(get_pos_in_lattice(i));
	set_pos_ab(x);
	set_in_LxLy(x);
	return unit_cell_index(x);
}

template<typename Type>
void System2D<Type>::set_pos_LxLy(Vector<double>& x) const {
	x=inv_LxLy_*x;
	set_pos(x);
}

template<typename Type>
void System2D<Type>::set_pos_ab(Vector<double>& x) const {
	x=inv_ab_*x;
	set_pos(x);
}
/*}*/

/*{private methods*/
template<typename Type>
void System2D<Type>::compute_TxTy(){
	//Tx_.set(this->n_,this->n_,0);
	//Ty_.set(this->n_,this->n_,0);
	//unsigned int tmp;
	//double t(1);
	//for(unsigned int j(0);j<Ly_;j++){
		//for(unsigned int i(0);i<Lx_-1;i++){
			//tmp = this->spuc_*(i + j*Lx_);
			//for(unsigned int k(0);k<this->spuc_;k++){
				//Tx_(tmp+k, tmp+k+this->spuc_) = t;
			//}
		//}
		//tmp = this->spuc_*((Lx_-1) + j*Lx_);
		//for(unsigned int k(0);k<this->spuc_;k++){
			//Tx_(tmp+k,this->spuc_*j*Lx_ + k) = this->bc_*t;
		//}
	//}
	//for(unsigned int i(0);i<Lx_;i++){
	//for(unsigned int j(0);j<Ly_-1;j++){
	//tmp = this->spuc_*(i + j*Lx_);
	//for(unsigned int k(0);k<this->spuc_;k++){
	//Ty_(tmp+k, tmp+this->spuc_*Lx_+k) = t;
	//}
	//}
	//tmp = this->spuc_*(i + (Ly_-1)*Lx_);
	//for(unsigned int k(0);k<this->spuc_;k++){
	//Ty_(tmp+k, this->spuc_*i+k) = this->bc_*t;
	//}
	//}
}

template<typename Type>
bool System2D<Type>::simple_diagonalization(){
	Vector<double> eval;
	Lapack<Type>(H_,false,(this->ref_(1)==1?'S':'H')).eigensystem(eval,true);
	for(unsigned int c(0);c<this->N_;c++){
		if(my::are_equal(eval(this->M_(c)),eval(this->M_(c)-1),1e-12)){
			std::cerr<<__PRETTY_FUNCTION__<<" : degenerate at the Fermi level"<<std::endl;
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
			if(my::are_equal(eval(i),eval(j),eq_prec_,eq_prec_)){
				std::cerr<<__PRETTY_FUNCTION__<<" : eigenvalue "<<i<<" and "<<j<<" degenerate"<<std::endl;
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

template<typename Type>
void System2D<Type>::set_pos(Vector<double>& x) const {
	double ip;
	x(0) = std::modf(x(0),&ip);
	x(1) = std::modf(x(1),&ip);
	if( my::are_equal(x(0),1,eq_prec_,eq_prec_) ){ x(0) = 0.0; }
	if( my::are_equal(x(1),1,eq_prec_,eq_prec_) ){ x(1) = 0.0; }
}

template<typename Type>
bool System2D<Type>::set_in_LxLy(Vector<double>& x) const {
	bool out_of_zone(false);
	double ip;
	x(0) = std::modf(x(0),&ip);
	if( x(0)<0 ){ x(0) += 1.0; out_of_zone = !out_of_zone; }
	if( my::are_equal(x(0),1,eq_prec_,eq_prec_) ){ x(0) = 0; out_of_zone = !out_of_zone; }
	if( ip>0 ){ out_of_zone = !out_of_zone; }

	x(1) = std::modf(x(1),&ip);
	if( x(1)<0 ){ x(1) += 1.0; out_of_zone = !out_of_zone; }
	if( my::are_equal(x(1),1,eq_prec_,eq_prec_) ){ x(1) = 0; out_of_zone = !out_of_zone; }
	if( ip>0 ){ out_of_zone = !out_of_zone ; }
	return out_of_zone;
}
/*}*/
#endif
