#include "PSpline.hpp"

PSpline::PSpline(unsigned int const& k):
	dim_(2),
	N_(0),
	phi_(&PSpline::phi1)
{
	switch(k){
		case 1: { phi_ = &PSpline::phi1; } break;
		case 2: { phi_ = &PSpline::phi2; } break;
		case 3: { phi_ = &PSpline::phi3; } break;
		case 4: { phi_ = &PSpline::phi4; } break;
		case 5: { phi_ = &PSpline::phi5; } break;
		case 6: { phi_ = &PSpline::phi6; } break;
		default: { std::cerr<<"PSpline::PSpline(unsigned int const& k) : unknown basis function"<<std::endl; } break;
	}
}

void PSpline::compute_weights(){
	Matrix<double> m(N_+dim_+1,N_+dim_+1);
	unsigned int k;
	for(unsigned int i(0);i<N_;i++){
		m(i,i) = 0.0;
		for(unsigned int j(i+1);j<N_;j++){
			m(i,j) = (this->*phi_)(sqrt((c_[i]-c_[j]).norm_squared()));
		}
		m(i,N_) = 1.0;
		k=0;
		for(unsigned int j(N_+1);j<N_+dim_+1;j++){
			m(i,j) = c_[i](k);
			k++;
		}
	}
	Lapack<double>(m,false,'S').inv();

	weights_.set(N_+dim_+1,0);
	for(unsigned int j(0);j<N_;j++){
		for(unsigned int i(0);i<N_;i++){
			weights_(i) += m(i,j)*y_[j];
		}
		for(unsigned int i(N_);i<N_+dim_+1;i++){
			weights_(i) += m(i,j)*y_[j];
		}
	}
}

double PSpline::extrapolate(Vector<double> const& x) const {
	double y(0);
	for(unsigned int i(0);i<N_;i++){
		y += weights_(i)*(this->*phi_)(sqrt((x-c_[i]).norm_squared()));
	}
	y += weights_(N_);
	for(unsigned int i(0);i<dim_;i++){ y += weights_(i+N_+1)*x(i); }
	return y;
}

void PSpline::add_data(Vector<double> const& c, double const& y){
	c_.push_back(c);
	y_.push_back(y);
	N_++;
}
