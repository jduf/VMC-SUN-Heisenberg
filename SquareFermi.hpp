#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

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
};


template<typename Type>
SquareFermi<Type>::SquareFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Square<Type>(1,1,1,"square-fermi")
{
	if(this->status_==1){
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
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(ll);
	double ey(ll);
	std::string color;

	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+this->filename_);
	unsigned int s;
	for(unsigned int i(0);i<this->Lx_;i++) {
		for(unsigned int j(0);j<this->Ly_;j++) {
			s = this->spuc_*(i+j*this->Lx_);
			x0 = i*ex;
			y0 = j*ey;
			nb = this->get_neighbourg(s);
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
		}
	}

	ps.frame(-0.5,-0.5,this->Lx_*ex-0.5,this->Ly_*ey-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,ex-0.5,ey-0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

template<typename Type>
void SquareFermi<Type>::check(){
	//unsigned int c(0);
	//unsigned int a(this->M_(c)-1);
	//unsigned int b(this->M_(c)-1);
	//Vector<double> eval;
	//do{b++;} while (b+1<this->n_ && are_equal(eval(b),eval(b-1)));
	//if(b!=this->M_(c)){ while(a>0 && are_equal(eval(a-1),eval(a))){a--;} }
	//std::cout<<a<<" "<<b<<std::endl;
	this->compute_H();
	this->plot_band_structure();
	this->degenerate_ = true;
}
/*}*/
#endif
