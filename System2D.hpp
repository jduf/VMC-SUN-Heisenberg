#ifndef DEF_SYSTEM2D
#define DEF_SYSTEM2D

#include "GenericSystem.hpp"
#include "List.hpp"

/*{*//*!Abstract class that provides many methods to create 2D clusters

	   If the parent System doesn't have an obs_[0], this class will declare
	   and set all attributes that are required to create a 2D cluster (sites'
	   position, cluster's boundaries, translation vectors...). Otherwise,
	   none of these attributes are declared and the construction is very
	   efficient because everything that is required to compute the
	   wavefunction is already contained in obs_[0].

	   Once this class has been created, it can be used to create few
	   observables (bond energies, long range correlation...)

	   It can also be used to study/display the porperties of the wavefunction.
	   For instance it can compute the band structure by doing a full
	   diagonalization of H+Tx+Ty. It also provides methods that draw into a
	   PSTricks file.
	   *//*}*/
template<typename Type>
class System2D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System2D(Matrix<double> const& cluster_vertex, Matrix<double> const& ab, unsigned int const& spuc, unsigned int const& z, unsigned int const& ndir, std::string const& filename);
		/*!Destructor*/
		virtual ~System2D();

	protected:
		Matrix<double> const cluster_vertex_;	//!< vertices of the cluster
		Vector<double>* equivalent_vertex_;		//!< equivalent points by a cluster translation
		Vector<double>* boundary_vertex_;		//!< vertices of the boundary
		Vector<double>* dir_nn_;				//!< vectors towards nearest neighbours
		Vector<double>* x_;						//!< position of the sites

		/*!Set which observable will be measured*/
		void create_obs(unsigned int const& which_obs);
		/*!Create the bond energy observables*/
		virtual void bond_energy();
		/*!Create the long range correlation observables*/
		virtual void compute_long_range_correlation();
		/*!Create the color occupation observables*/
		virtual void color_occupation();

		/*!Check if the unit cell can correctly fit inside the cluster*/
		bool unit_cell_allowed();
		/*!Returns the index of the position x the unit cell basis (a,b)*/
		virtual unsigned int unit_cell_index(Vector<double> const& x) const = 0;

		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Returns the index of the site at position x (if outside the cluster, it is first reset inside)*/
		unsigned int site_index(Vector<double> x) const;
		/*!Returns true if x is outside the cluster*/
		bool pos_out_of_lattice(Vector<double> const& x) const;
		/*!Returns true if the segment (x0,x1) cross a boundary (not the cluster border)*/
		bool cross_boundary(Vector<double> const& x0, Vector<double> const& x1) const;

		/*!Plots the band structure E(px,py)*/
		void plot_band_structure();

		/*!Returns a matrix containing the vertices of the unit cell*/
		Matrix<double> draw_unit_cell(double const& xshift=0, double const& yshift=0) const;
		/*!Returns a matrix containing the vertices of the boundary*/
		Matrix<double> draw_boundary(bool const& full_boundary) const;
		/*!Draws the long range correlations contained int O in the PSTricks file*/
		void draw_long_range_correlation(PSTricks& ps, Observable const& O) const;
		/*!Computes and writes the flux per plaquette in the PSTricks file*/
		void draw_flux_per_plaquette(PSTricks& ps, unsigned int s0, Vector<double> x, double const& xd, double const& yd, unsigned int const& jj, unsigned int const& jp, unsigned int const& jn) const;

	private:
		Matrix<double> const ab_;//!< basis vectors of the unit cell  ((a_1,b_1),(a_2,b_2))
		Matrix<double> inv_ab_;	 //!< inverse of the matrix ab_
		Matrix<Type> Tx_;		 //!< translation operator along x-axis
		Matrix<Type> Ty_;		 //!< translation operator along y-axis
		Vector<double> px_;		 //!< eigenvalue of Tx
		Vector<double> py_;		 //!< eigenvalue of Ty
		Vector<double> e_;		 //!< eigenvalue of H_

		/*!Returns the index of the site i in the unit cell*/
		unsigned int site_index_to_unit_cell_index(unsigned int const& i) const;
		/*!Resets x so pos_out_of_lattice(x) returns true*/
		virtual bool reset_pos_in_lattice(Vector<double>& x) const = 0;
		/*!Get the vector that separates the site i from its neighbourg in the direction d*/
		virtual Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const = 0;

		/*!Computes the translation operators*/
		void compute_TxTy();
		/*!Diagonalizes H_+Tx_+Ty_ => compute the band structure E(px,py)*/
		bool full_diagonalization();
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

			if( this->ref_(4) && cluster_vertex_.ptr() ){
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
/*{observables*/
template<typename Type>
void System2D<Type>::create_obs(unsigned int const& which_obs){
	switch(which_obs){
		case 0: { for(unsigned int i(1);i<4;i++){ create_obs(i); } }break;
		case 1: { bond_energy(); }break;
		case 2: { compute_long_range_correlation(); }break;
		case 3: { color_occupation(); }break;
		default:{
					std::cerr<<__PRETTY_FUNCTION__<<" : unknown observable "<<which_obs<<std::endl;
					std::cerr<<"Available observables are :"<<std::endl;
					std::cerr<<" + Bond energy            : 1"<<std::endl;
					std::cerr<<" + Long range correlation : 2"<<std::endl;
					std::cerr<<" + Color occupation       : 3"<<std::endl;
				}
	}
}

template<typename Type>
void System2D<Type>::bond_energy(){
	unsigned int idx(this->obs_.size());
	this->obs_.push_back(Observable("Bond energy",1,this->z_*this->spuc_/2,this->obs_[0].nlinks()));
	this->obs_[idx].remove_links();
}

template<typename Type>
void System2D<Type>::compute_long_range_correlation(){
	Vector<double>* dx(new Vector<double>[this->n_]);
	for(unsigned int i(0);i<this->n_;i++){ dx[i] = this->x_[i]-this->x_[0]; }

	this->obs_.push_back(Observable("Long range correlations",2,this->n_,this->n_*this->n_));
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			this->obs_.back()(i*this->n_+j,0) = i;
			this->obs_.back()(i*this->n_+j,1) = this->site_index(this->x_[i]+dx[j]);
			this->obs_.back()(i*this->n_+j,2) = j;
		}
	}
	delete[] dx;
}

template<typename Type>
void System2D<Type>::color_occupation(){
	this->obs_.push_back(Observable("Color occupation",3,this->N_*this->spuc_,Matrix<int>(this->n_,this->N_),this->n_/this->spuc_));
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->N_;j++){
			this->obs_.back()(i,j) = this->obs_[0](2*i,5)*this->N_+j;
		}
	}
}
/*}*/

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
	alpha = std::abs(std::modf(inv_ab_(0,0)*(equivalent_vertex_[1](0)-equivalent_vertex_[0](0))+inv_ab_(0,1)*(equivalent_vertex_[1](1)-equivalent_vertex_[0](1)),&ip));
	beta  = std::abs(std::modf(inv_ab_(1,0)*(equivalent_vertex_[1](0)-equivalent_vertex_[0](0))+inv_ab_(1,1)*(equivalent_vertex_[1](1)-equivalent_vertex_[0](1)),&ip));
	if( (!my::are_equal(alpha,0.0) && !my::are_equal(alpha,1.0)) || (!my::are_equal(beta,0.0) && !my::are_equal(beta,1.0)) ){
		std::cerr<<__PRETTY_FUNCTION__<<" : unit cell doesn't fit into the cluster (not sure) alpha="<<alpha<<" beta="<<beta<<std::endl;
		return false;
	}

	alpha = std::abs(std::modf(inv_ab_(0,0)*(equivalent_vertex_[2](0)-equivalent_vertex_[0](0))+inv_ab_(0,1)*(equivalent_vertex_[2](1)-equivalent_vertex_[0](1)),&ip));
	beta  = std::abs(std::modf(inv_ab_(1,0)*(equivalent_vertex_[2](0)-equivalent_vertex_[0](0))+inv_ab_(1,1)*(equivalent_vertex_[2](1)-equivalent_vertex_[0](1)),&ip));
	if( (!my::are_equal(alpha,0.0) && !my::are_equal(alpha,1.0)) || (!my::are_equal(beta,0.0) && !my::are_equal(beta,1.0)) ){
		std::cerr<<__PRETTY_FUNCTION__<<" : unit cell doesn't fit into the cluster (not sure) alpha="<<alpha<<" beta="<<beta<<std::endl;
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
		nb(d,2) = cross_boundary(x_[i],nn[d]);
		reset_pos_in_lattice(nn[d]);
	}

	unsigned int j(0);
	do {
		for(unsigned int d(0);d<dir.size();d++){
			if(my::are_equal(x_[j],nn[dir[d]],this->eq_prec_,this->eq_prec_)){
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
unsigned int System2D<Type>::site_index(Vector<double> x) const {
	reset_pos_in_lattice(x);
	for(unsigned int i(0);i<this->n_;i++){
		if(my::are_equal(x_[i],x,this->eq_prec_,this->eq_prec_)){ return i; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unkown site at position "<<x<<std::endl;
	return this->n_;
}

template<typename Type>
bool System2D<Type>::pos_out_of_lattice(Vector<double> const& x) const {
	return !my::in_polygon(cluster_vertex_.row(),cluster_vertex_.ptr(),cluster_vertex_.ptr()+cluster_vertex_.row(),x(0),x(1));
}

template<typename Type>
bool System2D<Type>::cross_boundary(Vector<double> const& x0, Vector<double> const& x1) const {
	bool does(false);
	if(my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[0].ptr(),boundary_vertex_[1].ptr()) || my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[2].ptr(),boundary_vertex_[3].ptr()) ){ does=!does; }
	if(my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[4].ptr(),boundary_vertex_[5].ptr()) || my::intersect(x0.ptr(),x1.ptr(),boundary_vertex_[6].ptr(),boundary_vertex_[7].ptr()) ){ does=!does; }
	return does;
}

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
			if(!my::are_equal(x,l.get()(0),this->eq_prec_,this->eq_prec_)){
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

/*{to draw the lattice*/
template<typename Type>
Matrix<double> System2D<Type>::draw_unit_cell(double const& xshift, double const& yshift) const {
	Matrix<double> tmp(4,2);
	tmp(0,0) = xshift;//  - (ab_(0,0)+ab_(0,1))/2.0;
	tmp(0,1) = yshift;//  - (ab_(1,0)+ab_(1,1))/2.0;
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

template<typename Type>
void System2D<Type>::draw_long_range_correlation(PSTricks& ps, Observable const& O) const {
	std::string color;
	double corr;
	double rescale(std::abs(0.25/O[1].get_x()));
	ps.cross(x_[0],0.25,"linecolor=black");
	ps.circle(x_[0],0.25,"linecolor=black");
	for(unsigned int i(1);i<this->n_;i++){
		corr = O[i].get_x();
		if(corr>0){ color = "blue"; }
		else      { color = "red"; }
		corr = sqrt(std::abs(corr*rescale));
		if(corr>1e-4 && std::abs(O[i].get_x())>O[i].get_dx()){ ps.circle(x_[i],corr,"fillstyle=solid,fillcolor="+color+",linecolor="+color); }
		else{ ps.cross(x_[i],0.1,"linecolor=black"); }
	}
}

template<typename Type>
void System2D<Type>::draw_flux_per_plaquette(PSTricks& ps, unsigned int s0, Vector<double> x, double const& xd, double const& yd, unsigned int const& jj, unsigned int const& jp, unsigned int const& jn) const {
	unsigned int s1;
	unsigned int j(0);
	double flux(0.0);
	double sign;
	unsigned long long a;
	unsigned long long b;
	do {
		x += dir_nn_[jj*j+jp];
		s1 = this->site_index(x);
		flux+= std::arg(-this->H_(s0,s1));
		s0 = s1;
	} while (++j<jn);
	flux /= M_PI;
	if(!my::are_equal(flux,0.0,this->eq_prec_,this->eq_prec_)){
		if(my::to_fraction(flux,a,b,sign) && b!=1){ ps.put(xd,yd,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}"); }
		else if((unsigned int)(my::chop(flux))%2){  ps.put(xd,yd,"\\tiny{1}"); }
	}
}
/*}*/
/*}*/

/*{private methods*/
template<typename Type>
unsigned int System2D<Type>::site_index_to_unit_cell_index(unsigned int const& i) const {
	double ip;
	Vector<double> x(inv_ab_*(x_[i]-x_[0]));
	x(0) = std::modf(x(0),&ip);
	x(1) = std::modf(x(1),&ip);
	if( x(0)<0 ){ x(0)+= 1.0; }
	if( x(1)<0 ){ x(1)+= 1.0; }
	if( my::are_equal(x(0),1.0,this->eq_prec_,this->eq_prec_) ){ x(0) = 0.0; }
	if( my::are_equal(x(1),1.0,this->eq_prec_,this->eq_prec_) ){ x(1) = 0.0; }
	return unit_cell_index(x);
}

template<typename Type>
void System2D<Type>::compute_TxTy(){
	Vector<double> x;
	Tx_.set(this->n_,this->n_,0);
	Ty_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->n_;i++){
		x = x_[i];
		x(0)+= ab_(0,0);
		x(1)+= ab_(1,0);
		Tx_(i,site_index(x)) = (cross_boundary(x_[i],x)?this->bc_:1);

		x = x_[i];
		x(0)+= ab_(0,1);
		x(1)+= ab_(1,1);
		Ty_(i,site_index(x)) = (cross_boundary(x_[i],x)?this->bc_:1);
	}
	//std::cout<<this->H_*Tx_-Tx_*this->H_<<std::endl;
	//std::cout<<this->H_*Ty_-Ty_*this->H_<<std::endl;
	//std::cout<<Tx_<<std::endl;
	//std::cout<<Ty_<<std::endl;
}

template<typename Type>
bool System2D<Type>::full_diagonalization(){
	compute_TxTy();
	Matrix<Type> M(this->H_);
	M += Tx_*Type(3.0);
	M += Ty_*Type(7.0);
	Vector<std::complex<double> > eval;
	Lapack<Type>(M,false,'G').eigensystem(eval,&this->evec_);

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(my::are_equal(eval(i),eval(j),this->eq_prec_,this->eq_prec_)){
				std::cerr<<__PRETTY_FUNCTION__<<" : eigenvalue "<<i<<" and "<<j<<" degenerate"<<std::endl;
				//return false;
			}
		}
	}
	Vector<unsigned int> index;
	e_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		my::display_progress(i,this->n_," compute E : ");
		e_(i) = this->projection(this->H_,i).real();
	}
	e_.sort(std::less_equal<double>(),index);

	Matrix<std::complex<double> > evec_tmp(this->evec_);
	Vector<std::complex<double> > eval_tmp(eval);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			std::swap(this->evec_(i,j),evec_tmp(i,index(j)));
		}
		std::swap(eval(i),eval_tmp(index(i)));
	}

	px_.set(this->n_);
	py_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		my::display_progress(i,this->n_," compute (px,py) : ");
		px_(i) = log(this->projection(Tx_,i)).imag();
		py_(i) = log(this->projection(Ty_,i)).imag();
	}
	return true;
}
/*}*/
#endif
