#include "SquareMu.hpp"

SquareMu::SquareMu(Parseur& P):
	Square<double>(P),
	mu_(P.get<double>("mu"))
{
	if(!P.status()){
		if(n_==Ly_*Lx_){
			if(study_system_){
				std::cerr<<"cs : SquareMu : will not create jdbin file"<<std::endl;
				compute_T(0);
				unsigned int color(P.get<bool>("color"));
				compute_band_structure(color);
				diagonalize_T('S');
				study_system();
			} else {
				for(unsigned int color(0);color<N_;color++){
					compute_T(color);
					diagonalize_T('S');
					for(unsigned int i(0);i<n_;i++){
						for(unsigned int j(0);j<m_;j++){
							EVec_(i+color*n_,j) = T_(i,j);
						}
					}
				}
				if(successful_){
					std::string filename("square-stripe");
					filename += "-N" + tostring(N_);
					filename += "-S" + tostring(n_);
					filename += "-" + tostring(Lx_) + "x" + tostring(Ly_);
					if(bc_ == 1){ filename += "-P";} 
					else { filename += "-A";}
					filename += "-mu" + tostring(mu_);
					save(filename);
				} else {
					std::cerr<<"SquareMu : degeneate"<<std::endl;
				}
			}
		} else {
			std::cerr<<"SquareMu : the cluster is not a square"<<std::endl;
		}
	}
}

SquareMu::~SquareMu(){}

void SquareMu::compute_T(unsigned int color){
	T_.set(n_,n_,0.0);
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*chemical potential*/
		if( (i-color) % N_ == 0 && i >= color){ T_(i,i) = mu_/2; }
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;  color++; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	T_ += T_.transpose();
}

void SquareMu::compute_P(){
	Px_.set(n_,n_,0.0);
	Py_.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - N_ ){Px_(i,i+N_) = 1; }
		else{ Px_(i,i-Ly_+N_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){
			if( (i+1) % Lx_ ){Py_(i,i+Lx_+1) = 1; }
			else { Py_(i,i+1) = bc_;}
		} else {
			if( (i+1) % Lx_ ) {  Py_(i,i-(Ly_-1)*Lx_+1) = bc_;}
			else { Py_(i,0) = bc_*bc_;}
		}
	}
}

void SquareMu::compute_band_structure(unsigned int color){
	compute_P();	

	//std::cout<<T_*Px_-Px_*T_<<std::endl;
	//std::cout<<T_*Py_-Py_*T_<<std::endl;

	Matrix<double> TP(T_+3.*Px_+7.*Py_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> kx(n_,1);
	Vector<double> ky(n_,1);
	Vector<double> E(n_,1);
	for(unsigned int i(0);i<n_;i++){
		kx(i) = log(projection(Px_,evec,i,i)).imag()/N_;
		ky(i) = log(projection(Py_,evec,i,i)).imag()-kx(i);
		E(i) = projection(T_,evec,i,i).real();
	}
	std::stringstream ss;
	ss<<color;
	std::string s("-"+ss.str());
	Gnuplot gp("spectrum" + s,"splot");
	gp.save_data("spectrum" + s,kx,ky,E);
	gp.add_plot_param(" ,\\\n");
	Vector<unsigned int> index(E.sort());
	gp.save_data("spectrum-sorted" + s,kx.sort(index).range(0,m_),ky.sort(index).range(0,m_),E.range(0,m_));
}

void SquareMu::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Stripe order : each color lives on its own sublattice");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("mu (chemical potential)",mu_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

void SquareMu::show(Matrix<double> const& T, unsigned int color){
	std::cout<<"T="<<std::endl;
	std::cout<<T<<std::endl;
	std::cout<<"favored sites :"<<std::endl;
	for(unsigned int i(0);i<Ly_;i++){
		for(unsigned int j(0);j<Lx_;j++){
			if(T(i+j*Ly_,i+j*Ly_)!=0){
				std::cout<<color<<" ";
			} else { 
				std::cout<<0<<" ";
			}
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}


