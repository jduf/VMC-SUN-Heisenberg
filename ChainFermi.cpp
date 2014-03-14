#include"ChainFermi.hpp"

ChainFermi::ChainFermi(unsigned int N, unsigned int n, unsigned int m, int bc):
	Chain<double>(N,n,m,bc,"chain-fermi")
{
	rst_.text("Spin ChainFermi, all the hopping parameters are real");
}

ChainFermi::~ChainFermi(){}

unsigned int ChainFermi::create(double x){
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
	if(degenerate_){ return 0; }
	else { return 1; }/*1st step successful*/
}

void ChainFermi::compute_T(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i< n_; i++){
		nb = get_neighbourg(i);
		T_(i,nb(0,0)) = nb(0,1)*t;
	}
	T_ += T_.transpose();
}

void ChainFermi::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainFermi::check(){
	bc_ = -1;
	compute_T();
	std::cout<<T_<<std::endl;
}

void ChainFermi::study(double E, double DeltaE, Vector<double> const& corr, std::string save_in){
	PSTricks ps(save_in,filename_,false);
	ps.add("\\begin{pspicture}(-1,-1)(13,7.75)%"+filename_);
	ps.put(10,7,"$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ $E\\pm\\Delta E="+tostring(E)+"\\pm"+tostring(DeltaE)+"$");
	std::string color;
	double ll(2.0); //link length
	unsigned int i(0);
	unsigned int j(0);
	for(unsigned int k(0);k<links_.row();k++){
		if(corr(k)>0){color="blue";}
		else{color="red";}
		ps.line("-",ll*i,j,ll*(i+1),j,"linewidth="+ tostring(corr(k))+"pt,linecolor="+color);
		ps.put(ll*i,j+0.2,tostring(k));
		ps.put(ll*(i+0.5),j-0.2,tostring(corr(k)));
		i++;
		if(i%int(sqrt(n_))==0){
			j++;
			i=0;
		}
	}
	ps.add("\\end{pspicture}");
}
