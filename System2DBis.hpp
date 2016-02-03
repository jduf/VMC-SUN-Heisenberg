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

		void diagonalize(bool simple);

		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Reset x so that it belongs to the lattice (Lx,Ly)*/
		bool pos_out_of_lattice(Vector<double> const& x) const;
		/*!Returns the index of the site at position x*/
		unsigned int find_index(Vector<double> const& x) const;
		/*!If x1 is outside the cluster, resets x1 inside and returns bc*/
		bool handle_boundary(Vector<double> const& x0, Vector<double>& x1) const;

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
		/*!Returns the index of the site i in the unit cell*/
		unsigned int get_site_in_unit_cell(unsigned int const& i) const;
};

/*{constructors*/
template<typename Type>
System2DBis<Type>::System2DBis(Matrix<double> const& lattice_corners, Matrix<double> const& ab, unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	lattice_corners_(lattice_corners),
	boundary_(NULL),
	dir_nn_(NULL),
	x_(NULL),
	ab_(ab),
	eq_prec_(1e-12)
{
	if(this->spuc_){
		inv_ab_.set(2,2);
		inv_ab_(0,0) = ab_(1,1);
		inv_ab_(1,0) =-ab_(1,0);
		inv_ab_(0,1) =-ab_(0,1);
		inv_ab_(1,1) = ab_(0,0);
		inv_ab_/=(ab_(0,0)*ab_(1,1)-ab_(1,0)*ab_(0,1));
		this->status_ = 2;

		if(!this->obs_.size() || !this->obs_[0].nlinks()){
			if(lattice_corners_.ptr()){
				boundary_ = new Vector<double>[4];
				dir_nn_ = new Vector<double>[this->z_];
				x_ = new Vector<double>[this->n_];		
				for(unsigned int i(0);i<this->n_;i++){ x_[i].set(2); }
				for(unsigned int i(0);i<this->z_;i++){ dir_nn_[i].set(2); }

				/*!check if the distance between each corners matches the
				 * vectors of the unit cell*/
				double alpha;
				double beta;
				double ip;
				bool fit(true);
				for(unsigned int i(0);i<lattice_corners_.row()-1;i++){
					alpha = std::modf(inv_ab_(0,0)*(lattice_corners_(i+1,0)-lattice_corners_(i,0))+inv_ab_(0,1)*(lattice_corners_(i+1,1)-lattice_corners_(i,1)),&ip);
					beta  = std::modf(inv_ab_(1,0)*(lattice_corners_(i+1,0)-lattice_corners_(i,0))+inv_ab_(1,1)*(lattice_corners_(i+1,1)-lattice_corners_(i,1)),&ip);
					if( !my::are_equal(alpha,0.0) ||!my::are_equal(beta,0.0) ){ fit = false; }
				}

				if(!fit){
					std::cerr<<__PRETTY_FUNCTION__<<" : it seems that the unit cell doesn't fit into the cluster (not sure)"<<std::endl;
					this->status_ = 3;
				}
			} else { std::cerr<<__PRETTY_FUNCTION__<<" : undefined geometry"<<std::endl; }
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the unit cell contains 0 site"<<std::endl; }
}

template<typename Type>
System2DBis<Type>::~System2DBis(){
	if(boundary_){ delete[] boundary_; }
	if(x_){ delete[] x_; }
	if(dir_nn_){ delete[] dir_nn_; }
}
/*}*/

/*{protected methods*/
template<typename Type>
void System2DBis<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_ = 1; } }
	else { if(full_diagonalization()){ this->status_ = 1; } }
}

template<typename Type>
void System2DBis<Type>::plot_band_structure(){
	if(full_diagonalization()){
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

		double min_e(e_.min());
		double max_e(e_.max());
		IOFiles bs(this->filename_+"-band-structure.dat",true);
		for(unsigned int i(0);i<this->n_;i++){
			a = std::make_shared<Vector<double> >(2+this->spuc_,666);
			b = NULL;
			(*a)(0) = my::chop(my::are_equal(std::abs(px_(i)),M_PI,1e-12)?-M_PI:px_(i));
			(*a)(1) = my::chop(my::are_equal(std::abs(py_(i)),M_PI,1e-12)?-M_PI:py_(i));

			bs<<(*a)(0)<<" "<<(*a)(1)<<" "<<e_(i)<<IOFiles::endl;

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

		IOFiles bsf(this->filename_+"-band-structure-formated.dat",true);
		l.set_target();
		double x(666);
		while(l.target_next()){
			if(!my::are_equal(x,l.get()(0))){
				x = l.get()(0);
				bsf<<IOFiles::endl;
			}
			bsf<<l.get()<<IOFiles::endl;
		}

		Gnuplot gp("./",this->filename_+"-band-structure");
		gp.range("x","-pi","pi");
		gp.range("y","-pi","pi");
		gp.range("z",min_e,max_e);
		gp += "set ticslevel 0";
		for(unsigned int i(0);i<this->spuc_;i++){
			gp+=std::string(!i?"splot":"     ")+" '"+this->filename_+"-band-structure-formated.dat' u 1:2:"+my::tostring(i+3)+" w l notitle"+(i+1==this->spuc_?"":",\\");
		}
		gp.save_file();
		gp.create_image(true,true);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : band structure not plotted"<<std::endl; }
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
		nb(d,1) = handle_boundary(x_[i],nn[d]);
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
	//std::cerr<<"-----"<<std::endl;
	//std::cerr<<i<<std::endl;
	//for(unsigned int d(0);d<this->z_;d++){
	//std::cerr<< nn[d]<<std::endl;
	//}
	//std::cerr<<nb<<std::endl;
	//}
	assert(j<this->n_+1);
	delete[] nn;
	return nb;
}

template<typename Type>
unsigned int System2DBis<Type>::find_index(Vector<double> const& x) const {
	for(unsigned int i(0);i<this->n_;i++){
		if(my::are_equal(x_[i],x,eq_prec_,eq_prec_)){ return i; }
	}
	return 1;
}

template<typename Type>
unsigned int System2DBis<Type>::get_site_in_unit_cell(unsigned int const& i) const {
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

template<typename Type>
bool System2DBis<Type>::handle_boundary(Vector<double> const& x0, Vector<double>& x1) const {
	bool bc(false);
	if(pos_out_of_lattice(x1)){
		if(my::intersect(x0.ptr(),x1.ptr(),boundary_[0].ptr(),boundary_[1].ptr()) || my::intersect(x0.ptr(),x1.ptr(),boundary_[2].ptr(),boundary_[3].ptr()) ){ bc=!bc; }
		if(my::intersect(x0.ptr(),x1.ptr(),boundary_[1].ptr(),boundary_[2].ptr()) || my::intersect(x0.ptr(),x1.ptr(),boundary_[3].ptr(),boundary_[0].ptr()) ){ bc=!bc; }
		reset_pos_in_lattice(x1);
	}
	return bc;
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
	//std::cout<<Tx_<<std::endl;
	//std::cout<<Ty_<<std::endl;
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
	Lapack<Type>(M,false,'G').eigensystem(eval,&evec_);

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(my::are_equal(eval(i),eval(j),eq_prec_,eq_prec_)){
				std::cerr<<__PRETTY_FUNCTION__<<" : eigenvalue "<<i<<" and "<<j<<" degenerate"<<std::endl;
				//return false;
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
