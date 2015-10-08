#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

template<typename Type>
class SquareFermi: public Square<Type>{
	public:
		SquareFermi(System const& s);
		~SquareFermi() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice(std::string const& path, std::string const& filename);

		unsigned int match_pos_in_ab(Vector<double> const& x) const { (void)(x); return 0;};
};

template<typename Type>
SquareFermi<Type>::SquareFermi(System const& s):
	System(s),
	Square<Type>(1,1,0,"square-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
template<typename Type>
void SquareFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);
	Matrix<int> nb;
	for(unsigned int i(0); i < this->n_; i++){
		nb = this->get_neighbourg(i);
		for(unsigned int j(0);j<2;j++){
			this->H_(i,nb(j,0)) = nb(j,1);
		}
	}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void SquareFermi<Type>::lattice(std::string const& path, std::string const& filename){
	std::string color("black");
	std::string linestyle("solid");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,filename);
	ps.begin(-9,-10,16,10,this->filename_);
	unsigned int s0;
	unsigned int s1;
	double t;
	for(unsigned int i(0);i<this->n_;i++) {
		s0 = this->link_types_[0](i,0);
		xy0 = this->get_pos_in_lattice(s0);
		this->set_pos_LxLy(xy0);
		xy0 = (this->LxLy_*xy0).chop();

		s1 = this->link_types_[0](i,1);
		xy1 = this->get_pos_in_lattice(s1);
		this->set_pos_LxLy(xy1);
		xy1 = (this->LxLy_*xy1).chop();

		if((xy0-xy1).norm_squared()<1.1){ linestyle = "solid"; }
		else {
			linestyle = "dashed";
			if(i%2 && xy1(1)<xy0(1)){
				xy1(0) = xy0(0);
				xy1(1) = xy0(1)+1.0;
			}
			if(!(i%2) && xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1.0;
				xy1(1) = xy0(1);
			}
			ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
		}

		t = this->H_(s0,s1);
		if(i%2){ ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }

		if(std::abs(t)>1e-4){
			if(t<0){ color = "red"; }
			else { color = "blue"; }

			xy0 = xy0.chop();
			xy1 = xy1.chop();
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);
		}
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

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=0;
	polygon(1,1)=1;
	polygon(2,0)=1;
	polygon(2,1)=1;
	polygon(3,0)=1;
	polygon(3,1)=0;
	ps.polygon(polygon,"linecolor=blue");

	ps.end(true,true,true);
}

template<typename Type>
void SquareFermi<Type>::check(){
	lattice("./","lattice");
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
