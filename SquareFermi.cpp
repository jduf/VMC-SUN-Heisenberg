#include "SquareFermi.hpp"

SquareFermi::SquareFermi(Container const& param):
	Square<double>(param,"square-fermi")
{}

SquareFermi::~SquareFermi(){}

void SquareFermi::compute_T(){
	double t(-1.0);
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0); i < n_; i++){
		neighbourg = get_neighbourg(i);
		for(unsigned int j(0);j<2;j++){
			T_(i,neighbourg(j)) = t;
		}
	}
	for(unsigned int i(0);i<BC_.row();i++){
		T_(BC_(i,0),BC_(i,1)) *= bc_;
	}
	T_ += T_.transpose();
}

void SquareFermi::create(double x){
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void SquareFermi::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("ref (wave function)",ref_);
	w("n (particles' number)",n_);
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("M (particles' number of each color)",M_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
}

void SquareFermi::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	Px.set(n_,n_,0.0);
	Py.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - 1 ){ Px(i,i+1) = 1; }
		else { Px(i,i+1-Lx_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){ Py(i,i+Lx_) = 1; }
		else { Py(i,i-(Ly_-1)*Lx_) = bc_; }
	}
}

void SquareFermi::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0);i<n_;i++){
		neighbourg = get_neighbourg(i);
		if((i+1) % Lx_ ){
			ps.line("-", i%Lx_, i/Ly_, neighbourg(0)%Lx_, neighbourg(0)/Ly_, "linewidth=1pt,linecolor=black");
		} else {
			ps.line("-", i%Lx_, i/Ly_, i%Lx_+1, neighbourg(0)/Ly_, "linewidth=1pt,linecolor=blue");
		}
		if( i+Lx_<this->n_){ 
			ps.line("-", i%Lx_, i/Ly_, neighbourg(1)%Lx_, neighbourg(1)/Ly_, "linewidth=1pt,linecolor=black");
		} else {
			ps.line("-", i%Lx_, i/Ly_, neighbourg(1)%Lx_, i/Ly_+1, "linewidth=1pt,linecolor=blue");
		}
	}

	double r(0.2);
	Vector<double> tmp(2);
	for(unsigned int i(0);i<n_;i++){
		ps.put(i%Lx_+r*1.2, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,0.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}

void SquareFermi::study(){
	compute_T();
	Matrix<double> Px;
	Matrix<double> Py;
	compute_P(Px,Py);
	//band_structure();
	lattice();
}
