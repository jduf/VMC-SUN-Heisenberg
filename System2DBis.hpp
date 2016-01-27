#ifndef DEF_SYSTEM2DBIS
#define DEF_SYSTEM2DBIS

#include "GenericSystem.hpp"
#include "List.hpp"

/*{Description*/
/*!The rules that the arguments of the constructor must obey are the folloing :
 *
 * This convention allows the use of local basis in the definitions of :
 *
 * + virtual Vector<double> get_pos_in_lattice(unsigned int const& i) const = 0;
 * + virtual unsigned int match_pos_in_ab(Vector<double> const& x) const = 0;
 */
/*}*/
template<typename Type>
class System2DBis: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System2DBis(Matrix<double> const& lattice_corners, Matrix<double> const& ab, unsigned int const& spuc, unsigned int const& z, std::string const& filename);
		/*!Destructor*/
		virtual ~System2DBis();

	protected:
		Matrix<Type> H_;					//!< matrix used to get the band structure
		Matrix<std::complex<double> > evec_;//!< eigenvectors of H+Tx+Ty
		Matrix<double> lattice_corners_;	//!< basis of the whole lattice
		Vector<double>* boundary_;
		Vector<double>* dir_nn_;
		Vector<double>* x_;
		Matrix<double> ab_;  		//!< basis of the unit cel
		double const eq_prec_;		//!< precision for equality (important for matchinf position in lattice)

		/*!Plots the band structure E(px,py)*/
		void plot_band_structure();
		/*!Creates the selection of optimal eigenvectors*/
		void select_eigenvectors();

		void diagonalize(bool simple);

		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Returns the index of the site i in the unit cell basis (a,b)*/
		unsigned int get_site_in_ab(unsigned int const& i) const;
		/*!Reset x so that it belongs to the lattice (Lx,Ly)*/
		bool pos_out_of_lattice(Vector<double> const& x) const;

	private:
		Matrix<Type> Tx_;		//!< translation operator along x-axis
		Matrix<Type> Ty_;		//!< translation operator along y-axis
		Vector<double> px_;		//!< eigenvalue of Tx
		Vector<double> py_;		//!< eigenvalue of Ty
		Vector<double> e_;		//!< eigenvalue of H_
		Matrix<double> inv_ab_;	//!< inverse of the matrix ab_

		/*!Computes the translation operators*/
		void compute_TxTy();
		/*!Diagonalizes H_*/
		bool simple_diagonalization();
		/*!Diagonalizes H_+T_ => compute the band structure E(p)*/
		bool full_diagonalization();
		/*!Evaluates the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);
		/*!Returns the index of the position x the unit cell basis (a,b)*/
		virtual unsigned int match_pos_in_ab(Vector<double> const& x) const = 0;
		/*!Resets x so it pos_out_of_lattice returns true*/
		virtual bool reset_pos_in_lattice(Vector<double>& x) const = 0;
		/*!Get the vector that separates the site i from its neighbourg in the direction d*/
		virtual Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const = 0;
};

/*{constructors*/
template<typename Type>
System2DBis<Type>::System2DBis(Matrix<double> const& lattice_corners, Matrix<double> const& ab, unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	lattice_corners_(lattice_corners),
	boundary_(new Vector<double>[4]),
	dir_nn_(new Vector<double>[this->z_]),
	x_(new Vector<double>[this->n_]),
	ab_(ab),
	eq_prec_(1e-12)
{
	if(lattice_corners_.size()){
		for(unsigned int i(0);i<this->n_;i++){ x_[i].set(2); }
		for(unsigned int i(0);i<this->z_;i++){ dir_nn_[i].set(2); }

		inv_ab_.set(2,2);
		inv_ab_(0,0) = ab_(1,1);
		inv_ab_(1,0) =-ab_(1,0);
		inv_ab_(0,1) =-ab_(0,1);
		inv_ab_(1,1) = ab_(0,0);
		inv_ab_/=(ab_(0,0)*ab_(1,1)-ab_(1,0)*ab_(0,1));

		if(this->spuc_){ this->status_--; }
		else { std::cerr<<__PRETTY_FUNCTION__<<" : the unit cell contains 0 site"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the number of site doesn't fit into the cluster"<<std::endl; }
}

template<typename Type>
System2DBis<Type>::~System2DBis(){
	if(dir_nn_){ delete[] dir_nn_; }
}
/*}*/

/*{protected methods*/
template<typename Type>
void System2DBis<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_--; } }
	else { if(full_diagonalization()){ this->status_--; } }
}

template<typename Type>
void System2DBis<Type>::plot_band_structure(){
	full_diagonalization();

	List<Vector<double> > l;
	std::shared_ptr<Vector<double> > a;
	List<Vector<double> >::Node* b;
	auto cmp = [](Vector<double> const& a, Vector<double> const& b){
		for(unsigned int i(0);i<2;i++){
			if(a(i) - b(i) > 0.0001){ return 0; }
			if(a(i) - b(i) <-0.0001){ return 1; }
		}
		return 2;
	};
	for(unsigned int i(0);i<this->n_;i++){
		a = std::make_shared<Vector<double> >(2+this->spuc_,666);
		b = NULL;
		(*a)(0) = my::chop(my::are_equal(std::abs(px_(i)),M_PI,1e-12)?-M_PI:px_(i));
		(*a)(1) = my::chop(my::are_equal(std::abs(py_(i)),M_PI,1e-12)?-M_PI:py_(i));
		if(l.find_in_sorted_list(a,b,cmp)){
			for(unsigned int j(2);j<2+this->spuc_;j++){
				if(e_(i)<(*b->get())(j)){
					std::swap(e_(i),(*b->get())(j));
				}
			}
		} else { 
			(*a)(2) = e_(i);
			l.set_target(b);
			l.add_after_target(a); 
		}
		l.set_target();
	}

	IOFiles spectrum("spectrum.dat",true);
	l.set_target();
	double x(666);
	while(l.target_next()){
		if(!my::are_equal(x,l.get()(0))){ 
			x = l.get()(0);
			spectrum<<IOFiles::endl;
		}
		spectrum<<l.get()<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.range("x","-pi","pi");
	gp.range("y","-pi","pi");
	gp.range("z","-5","5");
	for(unsigned int i(0);i<this->spuc_;i++){
		gp+=std::string(!i?"splot":"     ")+" 'spectrum.dat' u 1:2:"+my::tostring(i+3)+" w l notitle"+(i+1==this->spuc_?"":",\\");
	}
	gp.save_file();
	gp.create_image(true,true);
}

template<typename Type>
void System2DBis<Type>::select_eigenvectors(){
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
Matrix<int> System2DBis<Type>::get_neighbourg(unsigned int const& i) const {
	/*!nn* are the nearest neighbours */
	Vector<double>* nn(new Vector<double>[this->z_]);
	Matrix<int> nb(this->z_,2);
	std::vector<unsigned int> dir(this->z_);
	for(unsigned int d(0);d<this->z_;d++){
		dir[d]= d;
		nn[d] = x_[i];
		nn[d]+= get_relative_neighbourg_position(i,d);
		nb(d,1) = 0;
		if(my::intersect(x_[i].ptr(), nn[d].ptr(), boundary_[0].ptr(), boundary_[1].ptr()) || my::intersect(x_[i].ptr(), nn[d].ptr(), boundary_[2].ptr(), boundary_[3].ptr()) ){ nb(d,1) = !nb(d,1); }
		if(my::intersect(x_[i].ptr(), nn[d].ptr(), boundary_[1].ptr(), boundary_[2].ptr()) || my::intersect(x_[i].ptr(), nn[d].ptr(), boundary_[3].ptr(), boundary_[0].ptr()) ){ nb(d,1) = !nb(d,1); }
		reset_pos_in_lattice(nn[d]);
	}

	unsigned int j(0);
	do {
		for(unsigned int d(0);d<dir.size();d++){
			if(my::are_equal(x_[j],nn[dir[d]],eq_prec_,eq_prec_)){
				nb(dir[d],0) = j;
				dir.erase(dir.begin()+d);
			}
		}
	} while(dir.size() && ++j<this->n_+1);
	//if(j>=this->n_+1){
	//std::cout<<"-----"<<std::endl;
	//std::cout<<i<<std::endl;
	//for(unsigned int d(0);d<this->z_;d++){
	//std::cout<< nn[d]<<std::endl;
	//}
	//std::cout<<nb<<std::endl;
	//}
	assert(j<this->n_+1);
	delete[] nn;
	return nb;
}

template<typename Type>
unsigned int System2DBis<Type>::get_site_in_ab(unsigned int const& i) const {
	double ip;
	Vector<double> x(inv_ab_*(x_[i]-x_[0]));
	x(0) = std::modf(x(0),&ip);
	x(1) = std::modf(x(1),&ip);
	if( x(0)<0 ){ x(0)+= 1.0; }
	if( x(1)<0 ){ x(1)+= 1.0; }
	if( my::are_equal(x(0),1.0,eq_prec_,eq_prec_) ){ x(0) = 0.0; }
	if( my::are_equal(x(1),1.0,eq_prec_,eq_prec_) ){ x(1) = 0.0; }
	return match_pos_in_ab(x);
}

template<typename Type>
bool System2DBis<Type>::pos_out_of_lattice(Vector<double> const& x) const {
	return !my::in_polygon(lattice_corners_.row(),lattice_corners_.ptr(),lattice_corners_.ptr()+lattice_corners_.row(),x(0),x(1));
}
/*}*/

/*{private methods*/
template<typename Type>
void System2DBis<Type>::compute_TxTy(){
	Tx_.set(this->n_,this->n_,0);
	Vector<double> x;
	for(unsigned int i(0);i<this->n_;i++){
		x = x_[i];
		x(0)+= ab_(0,0);
		x(1)+= ab_(1,0);
		reset_pos_in_lattice(x);
		for(unsigned int j(0);j<this->n_;j++){
			if(my::are_equal(x,x_[j],eq_prec_,eq_prec_)){ Tx_(i,j) = 1; j=this->n_; }
		}
	}

	Ty_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->n_;i++){
		x = x_[i];
		x(0)+= ab_(0,1);
		x(1)+= ab_(1,1);
		reset_pos_in_lattice(x);
		for(unsigned int j(0);j<this->n_;j++){
			if(my::are_equal(x,x_[j],eq_prec_,eq_prec_)){ Ty_(i,j) = 1; j=this->n_; }
		}
	}
	//std::cout<<H_*Tx_-Tx_*H_<<std::endl;
	//std::cout<<H_*Ty_-Ty_*H_<<std::endl;
}

template<typename Type>
bool System2DBis<Type>::simple_diagonalization(){
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
bool System2DBis<Type>::full_diagonalization(){
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
std::complex<double> System2DBis<Type>::projection(Matrix<Type> const& O, unsigned int const& idx){
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
