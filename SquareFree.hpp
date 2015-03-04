#ifndef DEF_SQUAREFREE
#define DEF_SQUAREFREE

#include "Square.hpp"

template<typename Type>
class SquareFree: public Square<Type>{
	public:
		SquareFree(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t, Vector<double> const& mu);
		~SquareFree(){}

		void create();
		void check();

	protected:
		void compute_H(unsigned int const& c);
		void lattice();
		Vector<double> const t_;
		Vector<double> const mu_;
};

template<typename Type>
SquareFree<Type>::SquareFree(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t, Vector<double> const& mu):
	System(ref,N,m,n,M,bc),
	Square<Type>(2,2,1,"square-fermi"),
	t_(t),
	mu_(mu)
{
	//std::cout<<t_<<std::endl;
	//std::cout<<mu_<<std::endl;
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("Free : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
template<typename Type>
void SquareFree<Type>::create(){
	this->E_.set(50,5,false);
	this->corr_.set(this->links_.row(),50,5,false);

	for(unsigned int c(0);c<this->N_;c++){
		this->compute_H(c);
		this->diagonalize(true);
		if(this->status_==1){
			for(unsigned int i(0);i<this->n_;i++){
				for(unsigned int j(0);j<this->M_(c);j++){
					this->EVec_[c](i,j) = this->H_(i,j);
				}
			}
			if(c!=this->N_-1){ this->status_++;}
		}
	}
}

template<typename Type>
void SquareFree<Type>::compute_H(unsigned int const& c){
	this->H_.set(this->n_,this->n_,0);
	for(unsigned int i(0); i < this->links_.row(); i++){
		this->H_(this->links_(i,0),this->links_(i,1)) = t_(i%t_.size());
	}
	for(unsigned int i(0); i < this->n_; i++){
		if(i % this->N_ == c){ this->H_(i,i) = mu_(c)/2.0; }
	}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void SquareFree<Type>::lattice(){
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
void SquareFree<Type>::check(){
	//unsigned int c(0);
	//unsigned int a(this->M_(c)-1);
	//unsigned int b(this->M_(c)-1);
	//Vector<double> eval;
	//do{b++;} while (b+1<this->n_ && are_equal(eval(b),eval(b-1)));
	//if(b!=this->M_(c)){ while(a>0 && are_equal(eval(a-1),eval(a))){a--;} }
	//std::cout<<a<<" "<<b<<std::endl;
	std::cout<<t_<<std::endl;
	std::cout<<mu_<<std::endl;
	this->compute_H(1);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i);j<this->n_;j++){
			if(this->H_(i,j)!=0){std::cout<<i<<" "<<j<<" "<<this->H_(i,j)<<std::endl;}
		}
	}
	//this->plot_band_structure();
	this->status_++;
	std::cout<<this->H_<<std::endl;
}
/*}*/
#endif
