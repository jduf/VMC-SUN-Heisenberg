#include "SquareFermi.hpp"

SquareFermi::SquareFermi(unsigned int N, unsigned int n, unsigned int m):
	Square<double>(N,n,m,"square-fermi")
{
	rst_.text("Fermi : all colors experience the same Hamiltonian");
}

SquareFermi::~SquareFermi(){}

void SquareFermi::compute_T(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<2;j++){
			T_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	T_ += T_.transpose();
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
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		if((i+1) % Lx_ ){
			ps.line("-", i%Lx_, i/Ly_, nb(0,0)%Lx_, nb(0,0)/Ly_, "linewidth=1pt,linecolor=black");
		} else {
			ps.line("-", i%Lx_, i/Ly_, i%Lx_+1, nb(0,0)/Ly_, "linewidth=1pt,linecolor=blue");
		}
		if( i+Lx_<this->n_){ 
			ps.line("-", i%Lx_, i/Ly_, nb(1,0)%Lx_, nb(1,0)/Ly_, "linewidth=1pt,linecolor=black");
		} else {
			ps.line("-", i%Lx_, i/Ly_, nb(1,0)%Lx_, i/Ly_+1, "linewidth=1pt,linecolor=blue");
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

void SquareFermi::check(){
	compute_T();
	std::cout<<T_.chop(1e-6)<<std::endl;
}
