#include "SquareFreeReal.hpp"

SquareFreeReal::SquareFreeReal(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t, Vector<double> const& mu):
	System(ref,N,m,n,M,bc),
	Square<double>(2,2,1,"square-fermi"),
	t_(t),
	mu_(mu)
{
	//std::cout<<t_<<std::endl;
	//std::cout<<mu_<<std::endl;
	if(status_==2){
		init_fermionic();

		system_info_.text("FreeReal : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFreeReal::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	for(unsigned int c(0);c<N_;c++){
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
			if(c!=N_-1){ status_++;}
		}
	}
}

void SquareFreeReal::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);
	for(unsigned int i(0); i < links_.row(); i++){
		H_(links_(i,0),links_(i,1)) = t_(i%t_.size());
	}
	for(unsigned int i(0); i < n_; i++){
		if(i % N_ == c){ H_(i,i) = mu_(c)/2.0; }
	}
	H_ += H_.transpose();
}
/*}*/

/*{method needed for checking*/
void SquareFreeReal::lattice(){
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
	for(unsigned int i(0);i<Lx_;i++) {
		for(unsigned int j(0);j<Ly_;j++) {
			s = spuc_*(i+j*Lx_);
			x0 = i*ex;
			y0 = j*ey;
			nb = get_neighbourg(s);
			ps.put(x0-0.2,y0+0.2,my::tostring(s));
			x1 = x0+ll;
			y1 = y0;
			if(H_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*x-link*/ ps.line("-",x0,y0,x1,y1, "linewidth=1pt,linecolor="+color);
			x1 = x0;
			y1 = y0+ll;
			if(H_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*y-link*/ ps.line("-",x0,y0,x1,y1, "linewidth=1pt,linecolor="+color);
			x0 += ll;
		}
	}

	ps.frame(-0.5,-0.5,Lx_*ex-0.5,Ly_*ey-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,ex-0.5,ey-0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void SquareFreeReal::check(){
	//unsigned int c(0);
	//unsigned int a(M_(c)-1);
	//unsigned int b(M_(c)-1);
	//Vector<double> eval;
	//do{b++;} while (b+1<n_ && my::are_equal(eval(b),eval(b-1)));
	//if(b!=M_(c)){ while(a>0 && my::are_equal(eval(a-1),eval(a))){a--;} }
	//std::cout<<a<<" "<<b<<std::endl;
	std::cout<<t_<<std::endl;
	std::cout<<mu_<<std::endl;
	compute_H(1);
	for(unsigned int i(0);i<n_;i++){
		for(unsigned int j(i);j<n_;j++){
			if(H_(i,j)!=0){std::cout<<i<<" "<<j<<" "<<H_(i,j)<<std::endl;}
		}
	}
	//plot_band_structure();
	status_++;
	std::cout<<H_<<std::endl;
}
/*}*/
