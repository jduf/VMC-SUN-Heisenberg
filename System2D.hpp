#ifndef DEF_SYSTEM2D
#define DEF_SYSTEM2D

#include "GenericSystem.hpp"
#include "List.hpp"

/*{Description*/
/*!
 *
 * + virtual Vector<double> get_pos_in_lattice(unsigned int const& i) const = 0;
 * + virtual unsigned int unit_cell_index(Vector<double> const& x) const = 0;
 */
/*}*/
template<typename Type>
class System2D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System2D(Matrix<double> const& cluster_vertex, Matrix<double> const& ab, unsigned int const& spuc, unsigned int const& z, unsigned int const& ndir, std::string const& filename);
		/*!Destructor*/
		virtual ~System2D();

	protected:
		Matrix<Type> H_;						//!< matrix used to get the band structure
		Matrix<std::complex<double> > evec_;	//!< eigenvectors of H+Tx+Ty
		Matrix<double> const cluster_vertex_;	//!< vertices of the cluster
		Vector<double>* equivalent_vertex_;		//!< equivalent points by a cluster translation
		Vector<double>* boundary_vertex_;		//!< vertices of the boundary 
		Vector<double>* dir_nn_;				//!< vectors towards nearest neighbours
		Vector<double>* x_;						//!< position of the sites
		double const eq_prec_ = 1e-12;			//!< precision for equality (important for matching position in lattice)

		/*!Plots the band structure E(px,py)*/
		void plot_band_structure();
		/*!Diagonalize the trial Hamiltonian H_*/
		void diagonalize(bool simple);

		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Reset x so that it belongs to the lattice (Lx,Ly)*/
		bool pos_out_of_lattice(Vector<double> const& x) const;
		/*!If x1 is outside the cluster, resets x1 inside and returns bc*/
		bool handle_boundary(Vector<double> const& x0, Vector<double>& x1) const;

		/*!Returns the index of the site at position x*/
		unsigned int site_index(Vector<double> const& x) const;
		/*!Returns the index of the position x the unit cell basis (a,b)*/
		virtual unsigned int unit_cell_index(Vector<double> const& x) const = 0;

		/*!Check if the unit cell can correctly fit inside the cluster*/
		bool unit_cell_allowed();
		/*!Returns a matrix containing the vertices of the unit cell*/
		Matrix<double> draw_unit_cell(double const& xshift=0, double const& yshift=0) const;
		/*!Returns a matrix containing the vertices of the boundary*/
		Matrix<double> draw_boundary(bool const& full_boundary) const;

	private:
		Matrix<double> const ab_;//!< basis vectors of the unit cell  ((a_1,b_1),(a_2,b_2))
		Matrix<double> inv_ab_;	 //!< inverse of the matrix ab_
		Matrix<Type> Tx_;		 //!< translation operator along x-axis
		Matrix<Type> Ty_;		 //!< translation operator along y-axis
		Vector<double> px_;		 //!< eigenvalue of Tx
		Vector<double> py_;		 //!< eigenvalue of Ty
		Vector<double> e_;		 //!< eigenvalue of H_

		/*!Computes the translation operators*/
		void compute_TxTy();
		/*!Diagonalizes H_*/
		bool simple_diagonalization();
		/*!Diagonalizes H_+T_ => compute the band structure E(p)*/
		bool full_diagonalization();
		/*!Evaluates the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);
		/*!Resets x so it pos_out_of_lattice returns true*/
		virtual bool reset_pos_in_lattice(Vector<double>& x) const = 0;
		/*!Get the vector that separates the site i from its neighbourg in the direction d*/
		virtual Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const = 0;
		/*!Returns the index of the site i in the unit cell*/
		unsigned int get_site_in_unit_cell(unsigned int const& i) const;
};

/*{constructor*/
template<typename Type>
System2D<Type>::System2D(Matrix<double> const& cluster_vertex, Matrix<double> const& ab, unsigned int const& spuc, unsigned int const& z, unsigned int const& ndir, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	cluster_vertex_(cluster_vertex),
	equivalent_vertex_(NULL),
	boundary_vertex_(NULL),
	dir_nn_(NULL),
	x_(NULL),
	ab_(ab)
{
	if(this->status_==3){
		if(this->spuc_ && ab_.ptr()){
			inv_ab_.set(2,2);
			inv_ab_(0,0) = ab_(1,1);
			inv_ab_(1,0) =-ab_(1,0);
			inv_ab_(0,1) =-ab_(0,1);
			inv_ab_(1,1) = ab_(0,0);
			inv_ab_/=(ab_(0,0)*ab_(1,1)-ab_(1,0)*ab_(0,1));

			if( (!this->obs_.size() || !this->obs_[0].nlinks()) && cluster_vertex_.ptr()){
				equivalent_vertex_= new Vector<double>[3];
				boundary_vertex_ = new Vector<double>[8];
				dir_nn_ = new Vector<double>[ndir];
				x_ = new Vector<double>[this->n_];
				for(unsigned int i(0);i<this->n_;i++){ x_[i].set(2); }
				for(unsigned int i(0);i<ndir;i++){ dir_nn_[i].set(2); }
			}
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : the unit cell contains 0 site"<<std::endl; }
	}
}

template<typename Type>
System2D<Type>::~System2D(){
	if(equivalent_vertex_){ delete[] equivalent_vertex_; }
	if(boundary_vertex_){ delete[] boundary_vertex_; }
	if(x_){ delete[] x_; }
	if(dir_nn_){ delete[] dir_nn_; }
}
/*}*/

/*{protected methods*/
template<typename Type>
void System2D<Type>::plot_band_structure(){
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
					if(e_(i)<(*b->get())(j)){ std::swap(e_(i),(*b->get())(j)); }
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
			if(!my::are_equal(x,l.get()(0),eq_prec_,eq_prec_)){
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
void System2D<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_ = 1; } }
	else { if(full_diagonalization()){ this->status_ = 1; } }
}

template<typename Type>
bool System2D<Type>::unit_cell_allowed(){
	/*{*//*!Solves a system of equation, check if the solution belongs to N^2

		   v_i = ab_*Gamma =>  Gamma = ab_^(-1)*v_i

		   where Gamma=(alpha,beta) and v_i is the vector linking two
		   equivalent vertices. If alpha and beta are integer, then the unit
		   cell matches the cluster.
		   *//*}*/
	double ip;
	double alpha;
	double beta;
	alpha = std::modf(inv_ab_(0,0)*(equivalent_vertex_[1](0)-equivalent_vertex_[0](0))+inv_ab_(0,1)*(equivalent_vertex_[1](1)-equivalent_vertex_[0](1)),&ip);
	beta  = std::modf(inv_ab_(1,0)*(equivalent_vertex_[1](0)-equivalent_vertex_[0](0))+inv_ab_(1,1)*(equivalent_vertex_[1](1)-equivalent_vertex_[0](1)),&ip);
	if( !my::are_equal(alpha,0.0) || !my::are_equal(beta,0.0) ){
		std::cerr<<__PRETTY_FUNCTION__<<" : unit cell doesn't fit into the cluster (not sure)"<<std::endl; 
		return false; 
	}

	alpha = std::modf(inv_ab_(0,0)*(equivalent_vertex_[2](0)-equivalent_vertex_[0](0))+inv_ab_(0,1)*(equivalent_vertex_[2](1)-equivalent_vertex_[0](1)),&ip);
	beta  = std::modf(inv_ab_(1,0)*(equivalent_vertex_[2](0)-equivalent_vertex_[0](0))+inv_ab_(1,1)*(equivalent_vertex_[2](1)-equivalent_vertex_[0](1)),&ip);
	if( !my::are_equal(alpha,0.0) || !my::are_equal(beta,0.0) ){
		std::cerr<<__PRETTY_FUNCTION__<<" : unit cell doesn't fit into the cluster (not sure)"<<std::endl; 
		return false; 
	}

	boundary_vertex_[0] = equivalent_vertex_[0]*2.0 - equivalent_vertex_[2];
	boundary_vertex_[1] = equivalent_vertex_[2]*2.0 - equivalent_vertex_[0];
	boundary_vertex_[2] = equivalent_vertex_[1]     + equivalent_vertex_[0] - equivalent_vertex_[2];
	boundary_vertex_[4] = equivalent_vertex_[0]*2.0 - equivalent_vertex_[1];
	boundary_vertex_[5] = equivalent_vertex_[1]*2.0 - equivalent_vertex_[0];
	boundary_vertex_[6] = equivalent_vertex_[2]     + equivalent_vertex_[0] - equivalent_vertex_[1];

	/*need to compute them like this so it works for square and honeycomb lattice*/
	boundary_vertex_[3] = boundary_vertex_[2] + (equivalent_vertex_[1]-boundary_vertex_[2])*3.0;
	boundary_vertex_[7] = boundary_vertex_[6] + (equivalent_vertex_[2]-boundary_vertex_[6])*3.0;

	return true;
}

template<typename Type>
Matrix<int> System2D<Type>::get_neighbourg(unsigned int const& i) const {
	/*!nn* are the nearest neighbours */
	Vector<double>* nn(new Vector<double>[this->z_]);
	Matrix<int> nb(this->z_,3);
	std::vector<unsigned int> dir(this->z_);
	for(unsigned int d(0);d<this->z_;d++){
		dir[d]= d;
		nn[d] = x_[i];
		nn[d]+= get_relative_neighbourg_position(i,d,nb(d,1));
		nb(d,2) = handle_boundary(x_[i],nn[d]);
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
unsigned int System2D<Type>::site_index(Vector<double> const& x) const {
	for(unsigned int i(0);i<this->n_;i++){
		if(my::are_equal(x_[i],x,eq_prec_,eq_prec_)){ return i; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unkown site at position "<<x<<std::endl;
	return this->n_;
}

template<typename Type>
unsigned int System2D<Type>::get_site_in_unit_cell(unsigned int const& i) const {
	double ip;
	Vector<double> x(inv_ab_*(x_[i]-x_[0]));
	x(0) = std::modf(x(0),&ip);
	x(1) = std::modf(x(1),&ip);
	if( x(0)<0 ){ x(0)+= 1.0; }
	if( x(1)<0 ){ x(1)+= 1.0; }
	if( my::are_equal(x(0),1.0,eq_prec_,eq_prec_) ){ x(0) = 0.0; }
	if( my::are_equal(x(1),1.0,eq_prec_,eq_prec_) ){ x(1) = 0.0; }
	return unit_cell_index(x);
}

template<typename Type>
bool System2D<Type>::pos_out_of_lattice(Vector<double> const& x) const {
	return !my::in_polygon(cluster_vertex_.row(),cluster_vertex_.ptr(),cluster_vertex_.ptr()+cluster_vertex_.row(),x(0),x(1));
}

template<typename Type>
bool System2D<Type>::handle_boundary(Vector<double> const& x0, Vector<double>& x1) const {
	bool bc(false);
	if(my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[0].ptr(),boundary_vertex_[1].ptr()) || my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[2].ptr(),boundary_vertex_[3].ptr()) ){ bc=!bc; }
	if(my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[4].ptr(),boundary_vertex_[5].ptr()) || my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[6].ptr(),boundary_vertex_[7].ptr()) ){ bc=!bc; }
	reset_pos_in_lattice(x1);
	return bc;
}

template<typename Type>
Matrix<double> System2D<Type>::draw_unit_cell(double const& xshift, double const& yshift) const {
	Matrix<double> tmp(4,2);
	tmp(0,0) = xshift   - (ab_(0,0)+ab_(0,1))/2.0;
	tmp(0,1) = yshift   - (ab_(1,0)+ab_(1,1))/2.0;
	tmp(1,0) = tmp(0,0) +  ab_(0,0);
	tmp(1,1) = tmp(0,1) +  ab_(1,0);
	tmp(2,0) = tmp(0,0) +  ab_(0,0)+ab_(0,1);
	tmp(2,1) = tmp(0,1) +  ab_(1,0)+ab_(1,1);
	tmp(3,0) = tmp(0,0) +  ab_(0,1);
	tmp(3,1) = tmp(0,1) +  ab_(1,1);
	return tmp;
}

template<typename Type>
Matrix<double> System2D<Type>::draw_boundary(bool const& full_boundary) const {
	Matrix<double> tmp(full_boundary?5:3,2);
	tmp(0,0)= equivalent_vertex_[1](0);
	tmp(0,1)= equivalent_vertex_[1](1);
	tmp(1,0)= equivalent_vertex_[0](0);
	tmp(1,1)= equivalent_vertex_[0](1);
	tmp(2,0)= equivalent_vertex_[2](0);
	tmp(2,1)= equivalent_vertex_[2](1);
	if(full_boundary){
		tmp(3,0)= equivalent_vertex_[2](0)+(equivalent_vertex_[1](0)-equivalent_vertex_[0](0));
		tmp(3,1)= equivalent_vertex_[2](1)+(equivalent_vertex_[1](1)-equivalent_vertex_[0](1));
		tmp(4,0)= equivalent_vertex_[1](0);
		tmp(4,1)= equivalent_vertex_[1](1);
	}
	return tmp;
}
/*}*/

/*{private methods*/
template<typename Type>
void System2D<Type>::compute_TxTy(){
	bool bc;
	Vector<double> x;
	Tx_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->n_;i++){
		x = x_[i];
		x(0)+= ab_(0,0);
		x(1)+= ab_(1,0);
		bc = handle_boundary(x_[i],x);
		for(unsigned int j(0);j<this->n_;j++){
			if(my::are_equal(x,x_[j],eq_prec_,eq_prec_)){ Tx_(i,j) = (bc?this->bc_:1); j=this->n_; }
		}
	}

	Ty_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->n_;i++){
		x = x_[i];
		x(0)+= ab_(0,1);
		x(1)+= ab_(1,1);
		bc = handle_boundary(x_[i],x);
		for(unsigned int j(0);j<this->n_;j++){
			if(my::are_equal(x,x_[j],eq_prec_,eq_prec_)){ Ty_(i,j) = (bc?this->bc_:1); j=this->n_; }
		}
	}
	//std::cout<<H_*Tx_-Tx_*H_<<std::endl;
	//std::cout<<H_*Ty_-Ty_*H_<<std::endl;
	//std::cout<<Tx_<<std::endl;
	//std::cout<<Ty_<<std::endl;
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
