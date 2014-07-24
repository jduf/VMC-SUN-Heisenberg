#include "SquareFermi.hpp"

SquareFermi::SquareFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Square<double>(1,1,1,"square-fermi")
{
	if(this->status_==1){
		init_fermionic();

		system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFermi::compute_H(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<2;j++){
			H_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	H_ += H_.transpose();
}

void SquareFermi::create(){
	compute_H();
	diagonalize_H(H_);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}
/*}*/

/*{method needed for checking*/
void SquareFermi::lattice(){
	PSTricks ps("./",filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		if((i+1) % Lx_ ){ ps.line("-", i%Lx_, i/Ly_, nb(0,0)%Lx_, nb(0,0)/Ly_, "linewidth=1pt,linecolor=black"); }
		else { ps.line("-", i%Lx_, i/Ly_, i%Lx_+1, nb(0,0)/Ly_, "linewidth=1pt,linecolor=blue"); }
		if( i+Lx_<this->n_){ ps.line("-", i%Lx_, i/Ly_, nb(1,0)%Lx_, nb(1,0)/Ly_, "linewidth=1pt,linecolor=black"); }
		else { ps.line("-", i%Lx_, i/Ly_, nb(1,0)%Lx_, i/Ly_+1, "linewidth=1pt,linecolor=blue"); }
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

void SquareFermi::check(){
	compute_H();
	std::cout<<H_.chop(1e-6)<<std::endl;
}
/*}*/
