#include"ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(unsigned int N, unsigned int n, unsigned int m, int bc):
	Chain<double>(N,n,m,bc,"chain-polymerized")
{
	rst_.nl();
	rst_.text("Spin chain, with different real hopping term.");
	rst_.text("For N colors and m particules per sites, every");
	rst_.text("N/m, there is a weaker bound, namely t-delta");
	rst_.text("instead of t+delta. (t=1,delta>0)");
}

ChainPolymerized::~ChainPolymerized(){}

unsigned int ChainPolymerized::create(double delta){
	filename_ += "-delta" + tostring(delta);
	delta_=delta;

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

void ChainPolymerized::compute_T(){
	/*!If t<0, delta<0 otherwise no polymerization occurs
	 * If t>0, delta>0 otherwise no polymerization occurs */
	double t(1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i += a_){
		for(unsigned int j(0); j<a_-1; j++){
			nb = get_neighbourg(i+j);
			T_(i+j,nb(0,0)) = t+delta_;
		}
		nb = get_neighbourg(i+a_-1);
		T_(i+a_-1,nb(0,0)) = nb(0,1)*(t-delta_);
	}
	T_ += T_.transpose();
}

void ChainPolymerized::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainPolymerized::save(Write& w) const{
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}

void ChainPolymerized::check(){
	delta_=0.1;
	compute_T();
	std::cout<<T_<<std::endl;
	std::cout<<links_<<std::endl;
}

void ChainPolymerized::study(double E, double DeltaE, Vector<double> const& corr, std::string save_in){
	PSTricks ps(save_in,filename_,true);
	ps.pdf();
	ps.add("\\begin{pspicture}(-1,-1)(13,7.75)%"+filename_);
	ps.put(-0.75,1.25,"$N="+tostring(N_)+"$","[bl]");
	ps.put(-0.75,0.75,"$m="+tostring(m_)+"$","[bl]");
	ps.put(-0.75,0.25,"$n="+tostring(n_)+"$","[bl]");
	ps.put(-0.75,-0.25,"$\\delta="+tostring(delta_)+"$","[bl]");
	ps.put(-0.75,-0.75,"$E\\pm\\Delta E="+tostring(E)+"\\pm"+tostring(DeltaE)+"$","[bl]");
	std::string color;
	double ll(11/sqrt(n_));//link length
	double x(1);
	double y(0.5);
	std::string style;
	std::string pt;
	for(unsigned int i(0);i<this->links_.row();i++){
		if(corr(i)>0){color="blue";}
		else{color="red";}
		pt = tostring(corr(i)); 
		if(std::abs(corr(i))<1e-7 ){ pt = "0";}
		if(std::abs(corr(i))>10 ){color="black"; pt="1";}
		ps.line("-",x,y,x+ll,y,"linewidth="+ pt +"pt,linecolor="+color);
		ps.put(x,y+0.2,"\\tiny{"+tostring(i)+"}");
		ps.put(x+0.5*ll,y-0.2,"\\tiny{"+pt+"}");
		x += ll;
		if((i+1)%int(sqrt(n_)) == 0){ 
			x = 1;
			y+= 7/sqrt(n_);
		}
	}
	ps.add("\\end{pspicture}");
}
