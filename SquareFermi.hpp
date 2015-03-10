#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"
#include "Rand.hpp"

template<typename Type>
class SquareFermi: public Square<Type>{
	public:
		SquareFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc);
		~SquareFermi(){}

		void create();
		void check();

	protected:
		void compute_H();
		void lattice();

		static Matrix<double> set_ab(){
			Matrix<double> tmp(2,2);
			tmp(0,0) = 1;
			tmp(1,0) = 0;
			tmp(0,1) = 0;
			tmp(1,1) = 1;
			return tmp;
		}
};

template<typename Type>
SquareFermi<Type>::SquareFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Square<Type>(SquareFermi<Type>::set_ab(),1,"square-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
template<typename Type>
void SquareFermi<Type>::compute_H(){
	double t(1.0);
	this->H_.set(this->n_,this->n_,0);
	Matrix<int> nb;
	for(unsigned int i(0); i < this->n_; i++){
		nb = this->get_neighbourg(i);
		for(unsigned int j(0);j<2;j++){
			this->H_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void SquareFermi<Type>::lattice(){
	Matrix<int> nb;
	std::string color("black");
	Vector<double> xy0_LxLy(2,0);
	Vector<double> xy1_LxLy(2,0);
	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-9,-10)(16,10)%"+this->filename_);
	for(unsigned int i(0);i<this->n_;i++) {
		xy0_LxLy = this->get_LxLy_pos(i);
		this->set_in_LxLy(xy0_LxLy);
		xy0_LxLy = (this->LxLy_*xy0_LxLy).chop();
		ps.put(xy0_LxLy(0)-0.20,xy0_LxLy(1)+0.15,tostring(i));
		nb = this->get_neighbourg(i);

		if(nb(0,1)<0){
			color = "red";
			xy1_LxLy = xy0_LxLy;
			xy1_LxLy(0) += 1.0;
			ps.put(xy1_LxLy(0)-0.20,xy1_LxLy(1)+0.15,tostring(nb(0,0)));
		} else {
			color = "black";
			xy1_LxLy = this->get_LxLy_pos(nb(0,0));
			this->set_in_LxLy(xy1_LxLy);
			xy1_LxLy = (this->LxLy_*xy1_LxLy).chop();
		}
		/*x-link*/ ps.line("-",xy0_LxLy(0),xy0_LxLy(1),xy1_LxLy(0),xy1_LxLy(1), "linewidth=1pt,linecolor="+color);

		if(nb(1,1)<0){
			color = "red";
			xy1_LxLy = xy0_LxLy;
			xy1_LxLy(1) += 1.0;
			ps.put(xy1_LxLy(0)-0.20,xy1_LxLy(1)+0.15,tostring(nb(1,0)));
		} else {
			color = "black";
			xy1_LxLy = this->get_LxLy_pos(nb(1,0));
			this->set_in_LxLy(xy1_LxLy);
			xy1_LxLy = (this->LxLy_*xy1_LxLy).chop();
		}
		/*y-link*/ ps.line("-",xy0_LxLy(0),xy0_LxLy(1),xy1_LxLy(0),xy1_LxLy(1), "linewidth=1pt,linecolor="+color);
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

template<typename Type>
void SquareFermi<Type>::check(){
	lattice();
	Matrix<int> nb;
	//Vector<int> dir(4,0);
	//Rand<unsigned int> rnd(0,3);
	//unsigned int d;
	//unsigned int p(0);
	//unsigned int iter(1000);
	//for(unsigned int i(0);i<iter;i++){
	//d=rnd.get();
	//nb=this->get_neighbourg(p);
	//p = nb(d,0);
	//dir(d)++;
	//}
	//std::cout<<p<<std::endl;
	//std::cout<<dir<<std::endl;
	//d=0;
	//p=0;
	//unsigned int i(0);
	//while(i<iter){
	//if(dir(d) != 0){
	//dir(d)--;
	//i++;
	//nb=this->get_neighbourg(p);
	//p = nb(d,0);
	//} else { d++; }
	//}
	//std::cout<<p<<std::endl;
	//
	//std::cout<<"######################"<<std::endl;
	//for(unsigned int i(0);i<this->n_;i++){
	//nb = this->get_neighbourg(i);
	//std::cout<<"i="<<i<<std::endl;
	//std::cout<<nb<<std::endl;
	//}
	std::cout<<"######################"<<std::endl;
	std::cout<<this->get_neighbourg(4)<<std::endl;

}
/*}*/
#endif
