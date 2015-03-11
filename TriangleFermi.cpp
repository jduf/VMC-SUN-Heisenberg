#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Triangle<double>(TriangleFermi::set_ab(),1,"triangle-fermi")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void TriangleFermi::compute_H(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<3;j++){
			H_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	H_ += H_.transpose();
}

void TriangleFermi::create(){
	compute_H();
	diagonalize(false);
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
void TriangleFermi::lattice(){
	Matrix<int> nb;
	std::string color("black");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-9,-10)(16,10)%"+this->filename_);
	for(unsigned int i(0);i<this->n_;i++) {
		xy0 = this->get_pos_in_lattice(i);
		this->set_pos_LxLy(xy0);
		this->set_in_basis(xy0);
		xy0 = (this->LxLy_*xy0).chop();
		ps.put(xy0(0)-0.20,xy0(1)+0.15,tostring(i));
		nb = this->get_neighbourg(i);

		if(nb(0,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,tostring(nb(0,0)));
		} else {
			color = "black";
			xy1 = this->get_pos_in_lattice(nb(0,0));
			this->set_pos_LxLy(xy1);
			this->set_in_basis(xy1);
			xy1 = (this->LxLy_*xy1).chop();
		}
		/*x-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color);

		if(nb(1,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(1) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,tostring(nb(1,0)));
		} else {
			color = "black";
			xy1 = this->get_pos_in_lattice(nb(1,0));
			this->set_pos_LxLy(xy1);
			this->set_in_basis(xy1);
			xy1 = (this->LxLy_*xy1).chop();
		}
		/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color);

		if(nb(2,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) -= 1.0;
			xy1(1) += 1.0;
			xy1 = xy1.chop();
			ps.put(xy1(0)-0.20,xy1(1)+0.15,tostring(nb(2,0)));
		} else {
			color = "black";
			xy1 = this->get_pos_in_lattice(nb(2,0));
			this->set_pos_LxLy(xy1);
			this->set_in_basis(xy1);
			xy1 = (this->LxLy_*xy1).chop();
		}
		/*xy-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color);
	}

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=this->LxLy_(0,0);
	polygon(1,1)=this->LxLy_(0,1);
	polygon(2,0)=this->LxLy_(0,0)+this->LxLy_(1,0);
	polygon(2,1)=this->LxLy_(0,1)+this->LxLy_(1,1);
	polygon(3,0)=this->LxLy_(1,0);
	polygon(3,1)=this->LxLy_(1,1);

	ps.polygon(polygon,"linecolor=green");
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void TriangleFermi::check(){
	//Matrix<int> nb;
	//std::cout<<"######################"<<std::endl;
	//for(unsigned int i(0);i<this->n_;i++){
		//nb = this->get_neighbourg(i);
		//std::cout<<"i="<<i<<std::endl;
		//std::cout<<nb<<std::endl;
	//}
	//std::cout<<"######################"<<std::endl;
	//nb = this->get_neighbourg(3);
	//std::cout<<nb<<std::endl;
	lattice();
}
/*}*/
