#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(Parseur& P):
	Triangle<double>(P,"triangle-fermi")
{
	if(!P.status()){
		if(P.get<bool>("study")){
			compute_T();
			//compute_P();
			//band_structure();
			lattice();
		} else {
			compute_T();
			diagonalize_T('S');
			for(unsigned int spin(0);spin<N_;spin++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<m_;j++){
						EVec_(i+spin*n_,j) = T_(i,j);
					}
				}
			}
			if(successful_){
				filename_ += "-N" + tostring(N_);
				filename_ += "-S" + tostring(n_);
				filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
				if(bc_ == 1){ filename_ += "-P";} 
				else { filename_ += "-A";}
				save();
			} else {
				std::cerr<<"TriangleFermi : degeneate"<<std::endl;
			}
		}
	}
}

TriangleFermi::~TriangleFermi(){}

void TriangleFermi::compute_T(){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
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
	T_ += T_.transpose();
}

void TriangleFermi::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

void TriangleFermi::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	Px.set(n_,n_,0.0);
	Py.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - 1 ){ Px(i,i+1) = 1; }
		else { Px(i,i+1-Lx_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){ Py(i,i+Lx_) = 1; }
		else { Py(i,i-(Ly_-1)*Lx_) = bc_; }
	}
}

void TriangleFermi::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	double e1(0.5);
	double e2(sqrt(3.0)/2.0);
	//double e1(0);/*will have exactly the same topology except for the link (0,N)*/
	//double e2(1);
	Vector<unsigned int> neighbourg;
	double x0, y0, x1, y1;
	std::string color;
	for(unsigned int i(0);i<n_;i++){
		neighbourg = get_neighbourg(i);
		x0=i%Lx_;
		y0=i/Ly_;

		color = "black";
		y1 = neighbourg(0)/Ly_;
		if((i+1) % Lx_ ){
			x1 = neighbourg(0)%Lx_;
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		if((i+1) % Lx_ ){
			x1 = neighbourg(1)%Lx_;
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		if( i+Lx_<this->n_){ 
			y1 = neighbourg(1)/Ly_;
		} else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		x1=neighbourg(2)%Lx_;
		if( i+Lx_<this->n_){ 
			y1 = neighbourg(2)/Ly_;
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
