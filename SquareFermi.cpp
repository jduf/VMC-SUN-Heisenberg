#include "SquareFermi.hpp"

SquareFermi::SquareFermi(Parseur& P):
	Square<double>(P)
{
	if(!P.status()){
		if(n_==Ly_*Lx_){
			if(study_system_){
				std::cerr<<"cs : SquareMu : will not create jdbin file"<<std::endl;
				compute_T();
				compute_band_structure();
				diagonalize_T('S');
				study_system();
			} else {
				compute_T();
				diagonalize_T('S');

				for(unsigned int spin(0);spin<N_;spin++){
					for(unsigned int i(0);i<n_;i++){
						for(unsigned int j(0);j<m_;j++){
							EVec_(i+spin*n_,j) = T_(i,j);
						}
					}
				}
				if(successful_){
					std::string filename("square-fermi");
					filename += "-N" + tostring(N_);
					filename += "-S" + tostring(n_);
					filename += "-" + tostring(Lx_) + "x" + tostring(Ly_);
					if(bc_ == 1){ filename += "-P";} 
					else { filename += "-A";}
					save(filename);
				} else {
					std::cerr<<"SquareFermi : degeneate"<<std::endl;
				}
			}
		} else {
			std::cerr<<"SquareFermi : the cluster is not a square"<<std::endl;
		}
	}
}

SquareFermi::~SquareFermi(){}

void SquareFermi::compute_T(){
	T_.set(n_,n_,0.0);
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
	}
	T_ += T_.transpose();
}

void SquareFermi::compute_P(){
	Px_.set(n_,n_,0.0);
	Py_.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - 1 ){ Px_(i,i+1) = 1; }
		else { Px_(i,i+1-Lx_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){ Py_(i,i+Lx_) = 1; }
		else { Py_(i,i-(Ly_-1)*Lx_) = bc_; }
	}
}

void SquareFermi::compute_band_structure(){
	compute_P();
	
	//std::cout<<T_*Px_-Px_*T_<<std::endl;
	//std::cout<<T_*Py_-Py_*T_<<std::endl;
	//std::cout<<Px_<<std::endl;
	
	Matrix<double> TP(T_+3.*Px_+7.*Py_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> kx(n_,1);
	Vector<double> ky(n_,1);
	Vector<double> E(n_,1);
	for(unsigned int i(0);i<n_;i++){
		kx(i) = log(projection(Px_,evec,i,i)).imag();
		ky(i) = log(projection(Py_,evec,i,i)).imag();
		E(i) = projection(T_,evec,i,i).real();
	}
	Gnuplot gp("spectrum","splot");
	gp.save_data("spectrum",kx,ky,E);
	gp.add_plot_param(" ,\\\n");
	Vector<unsigned int> index(E.sort());
	gp.save_data("spectrum-sorted",kx.sort(index).range(0,m_),ky.sort(index).range(0,m_),E.range(0,m_));
}

void SquareFermi::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}
