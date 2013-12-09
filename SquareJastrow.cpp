#include "SquareJastrow.hpp"

SquareJastrow::SquareJastrow(Parseur& P):
	Square<double>(P,"square-Jastrow"),
	nu_(P.get<double>("nu"))
{
	if(!P.status()){
		//if(P.get<bool>("study")){
			//compute_T();
			//band_structure();
		//} else {
			//compute_T();
			//diagonalize_T('S');
			//EVec_.set(N_*n_,n_,0.0);
			//for(unsigned int i(0);i<N_*n_;i++){
				//for(unsigned int j(0);j<n_;j++){
					//EVec_(i,j) = T_(i,j);
				//}
			//}
			//
			//if(successful_){
				//std::cout<<"prob"<<std::endl;
			//}
			////if(successful_){
				//filename_ += "-N" + tostring(N_);
				//filename_ += "-S" + tostring(n_);
				//filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
				//if(bc_ == 1){ filename_ += "-P";} 
				//else { filename_ += "-A";}
				//filename_ += "-AF+" + tostring(P.get<double>("AF"));
//
				//save();
			////} else {
				////std::cerr<<"SquareJastrow : degeneate"<<std::endl;
			////}
		//}

		filename_ += "-N" + tostring(N_);
		filename_ += "-S" + tostring(n_);
		filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		if(bc_ == 1){ filename_ += "-P";} 
		else { filename_ += "-A";}
		filename_ += "-nu+" + tostring(nu_);

		save();
	} else {
		std::cerr<<"SquareJastrow : need to provide nu"<<std::endl;
	}
}

SquareJastrow::~SquareJastrow(){}

void SquareJastrow::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Staggered magnetic field, Becca's idea to mimic an on site chemical potential");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("nu (jastrow coefficient)",nu_);
	w("sts (connected sites)",sts_);
}

//void SquareJastrow::compute_T(){
//double t(0);
//Matrix<double> tmpT(n_,n_,0.0);
//Matrix<double> tmpAF(n_,n_,0.0);
//T_.set(N_*n_,N_*n_,0.0);
//for(unsigned int i(0); i < n_; i++){
//tmpAF(i,i) = AF_; 
///*horizontal hopping*/
//if( (i+1) % Lx_ ){ tmpT(i,i+1) = t; AF_ *= -1;}
//else { tmpT(i+1-Lx_,i) = bc_*t;}
///*vertical hopping*/
//if( i+Lx_ < n_ ){  tmpT(i,i+Lx_) = t; } 
//else { tmpT(i-(Ly_-1)*Lx_,i) = bc_*t;}
//}
//
//for(unsigned int color(0); color<N_;color++){
//unsigned int c(color*n_);
//for(unsigned int i(0);i<n_;i++){
//for(unsigned int j(0);j<n_;j++){
//T_(c+i,c+j) = tmpT(i,j);
//}
//}
//}
//for(unsigned int color(1); color<N_;color++){
//unsigned int c(color*n_);
//for(unsigned int i(0);i<n_;i++){
//for(unsigned int j(0);j<n_;j++){
//T_(i,c+j) = tmpAF(i,j);
//}
//}
//}
//T_ += T_.transpose();
//}


//void SquareJastrow::compute_P(){
//Matrix<double> tmpPx(n_,n_,0.0);
//Matrix<double> tmpPy(n_,n_,0.0);
//Px_.set(N_*n_,N_*n_,0.0);
//Py_.set(N_*n_,N_*n_,0.0);
//for(unsigned int i(0); i < n_; i++){
///*horizontal hopping*/
//if( (i % Ly_)  < Ly_ - N_ ){tmpPx(i,i+N_) = 1; }
//else{ tmpPx(i,i-Ly_+N_) = bc_; }
///*vertical hopping*/
//if( i+Lx_ < n_ ){
//if( (i+1) % Lx_ ){tmpPy(i,i+Lx_+1) = 1; }
//else { tmpPy(i,i+1) = bc_;}
//} else {
//if( (i+1) % Lx_ ) { tmpPy(i,i-(Ly_-1)*Lx_+1) = bc_;}
//else { tmpPy(i,0) = bc_*bc_;}
//}
//}
//
//for(unsigned int color(0); color<N_;color++){
//unsigned int c(color*n_);
//for(unsigned int i(0);i<n_;i++){
//for(unsigned int j(0);j<n_;j++){
//Px_(c+i,c+j) = tmpPx(i,j);
//Py_(c+i,c+j) = tmpPy(i,j);
//}
//}
//}
//}

//void SquareJastrow::band_structure(){
//compute_P();
//
////std::cout<<T_*Px_-Px_*T_<<std::endl;
////std::cout<<T_*Py_-Py_*T_<<std::endl;
//
//Matrix<double> TP(T_+3.*Px_+7.*Py_);
//Vector<std::complex<double> > eval;
//Matrix<std::complex<double> > evec;
//Lapack<double> ES(&TP,false,'G');
//ES.eigensystem(&eval,&evec);
//Vector<double> kx(N_*n_);
//Vector<double> ky(N_*n_);
//Vector<double> E(N_*n_);
//for(unsigned int i(0);i<N_*n_;i++){
//kx(i) = log(projection(Px_,evec,i,i)).imag()/N_;
//ky(i) = log(projection(Py_,evec,i,i)).imag()-kx(i);
//E(i) = projection(T_,evec,i,i).real();
//}
//save_band_structure(kx,ky,E);
//}

//void SquareJastrow::lattice(){ }
