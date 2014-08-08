#include "SquareFermi.hpp"

SquareFermi::SquareFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Square<double>(1,1,1,"square-fermi")
{
	if(status_==1){
		init_fermionic();

		system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFermi::compute_H(){
	double t(1.0);
	H_.set(n_,n_,0);
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
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

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
	Matrix<int> nb;
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(ll);
	double ey(ll);
	std::string color;

	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+filename_);
	unsigned int s;
	for(unsigned int i(0);i<this->Lx_;i++) {
		for(unsigned int j(0);j<this->Ly_;j++) {
			s = this->spuc_*(i+j*this->Lx_);
			x0 = i*ex;
			y0 = j*ey;
			for(unsigned int k(0);k<spuc_;k++){
				nb = get_neighbourg(s);
				ps.put(x0-0.2,y0+0.2,tostring(s));
				x1 = x0+ll;
				y1 = y0;
				if(this->H_(s,nb(0,0))>0){ color = "green"; }
				else { color = "blue"; }
				/*x-link*/ ps.line("-",x0,y0,x1,y1, "linewidth=1pt,linecolor="+color);
				x1 = x0;
				y1 = y0+ll;
				if(this->H_(s,nb(1,0))>0){ color = "green"; }
				else { color = "blue"; }
				/*y-link*/ ps.line("-",x0,y0,x1,y1, "linewidth=1pt,linecolor="+color);
				x0 += ll;
				s++;
			}
		}
	}

	ps.frame(-0.5,-0.5,Lx_*ex-0.5,Ly_*ey-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,ex-0.5,ey-0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void SquareFermi::check(){
	Matrix<int> nb;
	for(unsigned int s(0);s<n_;s++){
		nb = get_neighbourg(s);
		for(unsigned int i(0);i<z_;i++){
			std::cout<<s<<" "<<nb(i,0)<<" "<<nb(i,1)<<std::endl;
		}
	}
	compute_H();
	lattice();
	std::cout<<H_.chop(1e-6)<<std::endl;
}
/*}*/
