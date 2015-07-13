#include "PSpline.hpp"

PSpline::PSpline(unsigned int const& k):
	basis_(k),
	dim_(2),
	N_(0),
	phi_(NULL)
{
	std::cerr<<"set dim"<<dim_<<std::endl;
}

void PSpline::compute_weights(){
	Vector<int> ipiv;
	Matrix<double> m(N_+dim_+1,N_+dim_+1);
	double rcn;
	while(!ipiv.ptr() && select_basis() ){
		unsigned int k;
//#pragma omp parallel for private(k)
		for(unsigned int i=0;i<N_;i++){
			m(i,i) = 0.0;
			for(unsigned int j(i+1);j<N_;j++){
				m(i,j) = (this->*phi_)(sqrt((c_[i]-c_[j]).norm_squared()));
			}
			k=0;
			for(unsigned int j(N_);j<N_+dim_;j++){
				m(i,j) = c_[i](k);
				k++;
			}
			m(i,N_+dim_) = 1.0;
		}
		for(unsigned int i(N_);i<N_+dim_+1;i++){
			for(unsigned int j(i);j<N_+dim_+1;j++){
				m(i,j) = 0.0;
			}
		}
		//Matrix<double> tmp(m);
		//std::cout<<m<<std::endl;
		//Lapack<double> inv_m(m,false,'S');
		//ipiv = inv_m.is_singular(rcn);
		//if(ipiv.ptr()){ 
			//inv_m.inv(ipiv);
			////for(unsigned int i(0);i<N_+dim_+1;i++){
				////for(unsigned int j(i);j<N_+dim_+1;j++){
					////tmp(j,i) = tmp(i,j);
				////}
			////}
			////std::cout<<(tmp*m).chop()<<std::endl;
		//} else {
			//std::cout<<"void PSpline::compute_weights() : can't invert the matrix with basis function "<<basis_<<" (rcn="<<rcn<<")"<<std::endl;
		//}
		Lapack<double>(m,false,'S').inv();
		ipiv.set(2);
	}

	weights_.set(N_+dim_+1,0);
	/*if I parallelize here, it bugs...*/
	//#pragma omp parallel for
	for(unsigned int i(0);i<N_+dim_+1;i++){
		for(unsigned int j(0);j<N_;j++){
			weights_(i) += m(i,j)*y_[j];
		}
	}
	for(unsigned int i(0);i<5;i++){
		std::cout<<y_[i]-extrapolate(c_[i])<<std::endl;
	}
}

double PSpline::extrapolate(Vector<double> const& x) const {
	double y(weights_.back());
	for(unsigned int i(0);i<N_;i++){
		y += weights_(i)*(this->*phi_)(sqrt((x-c_[i]).norm_squared()));
	}
	for(unsigned int i(0);i<dim_;i++){ y += weights_(i+N_)*x(i); }
	return y;
}

void PSpline::add_data(Vector<double> const& c, double const& y){
	c_.push_back(c);
	y_.push_back(y);
	N_++;
}

void PSpline::add_data(unsigned int const& i, Vector<double> const& c, double const& y){
	c_[i]=c;
	y_[i]=y;
}

void PSpline::set(unsigned int const& N){ 
	c_.clear();
	y_.clear();
	if(N){
		c_.resize(N);
		y_.resize(N);
		N_=N;
	} else { N_=0; }
}

bool PSpline::select_basis(){
	basis_++;
	std::cout<<"void PSpline::select_basis() : used basis function "<<basis_<<std::endl;
	switch(basis_){
		case 1: { phi_ = &PSpline::phi1; } break;
		case 2: { phi_ = &PSpline::phi2; } break;
		case 3: { phi_ = &PSpline::phi3; } break;
		case 4: { phi_ = &PSpline::phi4; } break;
		case 5: { phi_ = &PSpline::phi5; } break;
		case 6: { phi_ = &PSpline::phi6; } break;
		case 7: { phi_ = &PSpline::phi7; } break;
		default: {
					 basis_ = 0;
					 std::cerr<<"PSpline::PSpline(unsigned int const& k) : unknown basis function"<<std::endl;
					 return false;
				 } break;
	}
	return true;
}
