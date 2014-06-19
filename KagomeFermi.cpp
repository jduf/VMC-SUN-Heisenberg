#include "KagomeFermi.hpp"

KagomeFermi::KagomeFermi(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref):
	System(N,n,m,bc,ref),
	Kagome<double>(1,1,3,"kagome-fermi")
{
	rst_.text("KagomeFermi : All hopping term are identical, therefore the unit cell contains only 3 sites");
}

KagomeFermi::~KagomeFermi(){}

void KagomeFermi::compute_T(){
	double t(1.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<2;j++){
			T_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	T_ += T_.transpose();
}

void KagomeFermi::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	std::cerr<<"KagomeFermi::compute_P : undefined method"<<Px<<Py<<std::endl;
	Px.set(n_,n_,0);
	Py.set(n_,n_,0);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			for(unsigned int k(0);k<3;k++){
				Px(3*i + 3*j*Lx_ + k, 3*i + 3*j*Lx_ + k+3) = 1.0;
			}
		}
		for(unsigned int k(0);k<3;k++){
			Px(3*(Lx_-1) + 3*j*Lx_ +k,3*j*Lx_ + k) = bc_;
		}
	}
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_-1;j++){
			for(unsigned int k(0);k<3;k++){
				Py(3*i + 3*j*Lx_ + k, 3*i + 3*(j+1)*Lx_ + k) = 1.0;
			}
		}
		for(unsigned int k(0);k<3;k++){
			Py(3*i + 3*(Ly_-1)*Lx_ + k, 3*i + k) = bc_;
		}
	}
}

void KagomeFermi::create(double const& x, unsigned int const& type){
	if(type!=2){std::cerr<<"KagomeFermi::create(double x, unsigned int const& type) : type unknown"<<x<<std::endl;}
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void KagomeFermi::check(){
	compute_T();
	Matrix<double> Px;
	Matrix<double> Py;
	compute_P(Px,Py);
	BandStructure<double> bs(T_,Px,Py);

}
