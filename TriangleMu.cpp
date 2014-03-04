#include "TriangleMu.hpp"

TriangleMu::TriangleMu(unsigned int N, unsigned int n, unsigned int m):
	Triangle<double>(N,n,m,"triangle-mu")
{}

TriangleMu::~TriangleMu(){}

void TriangleMu::compute_T(unsigned int alpha){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*chemical potential*/
		/*if i(+)alpha stripe order perpenticular to the diagonal hopping*/
		/*if i(-)alpha stripe order parallel to the diagonal hopping*/
		if( (i+alpha) % N_ == 0 && i >= alpha){ T_(i,i) = mu_/2; }
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t; alpha++;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
		/*diagonal hopping*/
		if( (i+1) % Lx_ && i+Lx_ < n_ ){  T_(i,i+Lx_+1) = t; } 
		else {
			if(i+1 < n_ ){
				if( !((i+1) % Lx_) ){ T_(i,i+1) = bc_*t;}/*x jump across boundary*/
				if( i+Lx_ >= n_ ){  T_(i-Lx_*(Ly_-1)+1,i) = bc_*t; }/*y jump across boundary*/
			} else {
				T_(0,n_-1) = bc_*bc_*t;
			}
		}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	T_ += T_.transpose();
}

void TriangleMu::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	Px.set(n_,n_,0.0);
	Py.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - N_ ){Px(i,i+N_) = 1; }
		else{ Px(i,i-Ly_+N_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){
			if( i % Lx_ ){Py(i,i+Lx_-1) = 1; }
			else {Py(i,i+2*Lx_-1) = bc_;}
		} else {
			if( i % Lx_ ) { Py(i,i-(Ly_-1)*Lx_-1) = bc_;}
			else { Py(i,Lx_-1) = bc_*bc_;}
		}
	}
}

void TriangleMu::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	double e1(0.5);
	double e2(sqrt(3.0)/2.0);
	//double e1(0);/*will have exactly the same topology except for the link (0,N)*/
	//double e2(1);
	Matrix<int> nb;
	double x0, y0, x1, y1;
	std::string color;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		x0=i%Lx_;
		y0=i/Ly_;

		color = "black";
		y1 = nb(0,0)/Ly_;
		if((i+1) % Lx_ ){
			x1 = nb(0,0)%Lx_;
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		if((i+1) % Lx_ ){
			x1 = nb(1,0)%Lx_;
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		if( i+Lx_<this->n_){ 
			y1 = nb(1,0)/Ly_;
		} else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		x1=nb(2,0)%Lx_;
		if( i+Lx_<this->n_){ 
			y1 = nb(2,0)/Ly_;
		} else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);
	}

	diagonalize_T('S');

	double r(0.2);
	Vector<double> ada(n_,0);
	double max(occupation_number(ada));
	Vector<double> tmp(2);
	for(unsigned int i(0);i<n_;i++){
		x0 = i%Lx_;
		y0 = i/Ly_;
		tmp(0) = round(ada(i),7);
		tmp(1) = round((max-ada(i))/max,7);
		ps.add("\\rput("+tostring(x0-y0*e1)+","+tostring(y0*e2)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,white}");
		ps.add("}");
		ps.put(x0-y0*e1, y0*e2+r*1.2, "\\color{"+color+"}{\\tiny{"+tostring(i)+"}}");
	}

	Matrix<double> xy(4,2);
	xy(0,0) = -0.5*e1;
	xy(0,1) = -0.5;
	xy(1,0) = Ly_*(-e1)-0.5*e1;
	xy(1,1) = Ly_*e2-0.5;
	xy(2,0) = Lx_-0.5*e1 - Ly_*e1;
	xy(2,1) = Ly_*e2-0.5;;
	xy(3,0) = Lx_-0.5*e1;
	xy(3,1) = -0.5;

	ps.polygon(xy,"linecolor=red");
	ps.add("\\end{pspicture}");
}

void TriangleMu::study(){
	unsigned int alpha(1);
	compute_T(alpha);
	//compute_P();
	//for(unsigned int i(0);i<n_;i++){
	//kx(i) = log(projection(Px_,evec,i,i)).imag()/N_;
	//ky(i) = log(projection(Py_,evec,i,i)).imag()+kx(i);
	//E(i) = projection(T_,evec,i,i).real();
	//}
	//band_structure();
	lattice();
	std::cerr<<"TriangleMu : TriangleMu() : alpha must be smaller than N_"<<std::endl;
}

void TriangleMu::create(double mu){
	mu_ = mu;
	for(unsigned int alpha(0);alpha<N_;alpha++){
		compute_T(alpha);
		diagonalize_T('S');
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+alpha*n_,j) = T_(i,j);
			}
		}
		T_.set(n_,n_,0.0);
	}
}

void TriangleMu::save(Write& w) const{
	GenericSystem<double>::save(w);
	w("mu (chemical potential)",mu_);
}

void TriangleMu::check(){
	compute_T(1);
	std::cout<<T_.chop(1e-6)<<std::endl;
}
