#include "Interpolation.hpp"

Interpolation::Interpolation(unsigned int const& dim):
	dim_(dim),
	N_(0),
	support_(0),
	phi_(NULL)
{
	std::cerr<<"eps "<<support_<<std::endl;
}

void Interpolation::compute_weights(){
	Vector<int> ipiv;
	Matrix<double> m(N_+dim_+1,N_+dim_+1);
	double rcn;
	unsigned int k;
	for(unsigned int i=0;i<N_;i++){
		for(unsigned int j(i);j<N_;j++){
			m(i,j) = (this->*phi_)(sqrt((c_[i]-c_[j]).norm_squared())/support_);
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
	Lapack<double> inv_m(m,false,'S');
	ipiv = inv_m.is_singular(rcn);
	if(ipiv.ptr()){ 
		inv_m.inv(ipiv);
		weights_.set(N_+dim_+1,0);
		for(unsigned int i(0);i<N_+dim_+1;i++){
			for(unsigned int j(0);j<N_;j++){
				weights_(i) += m(i,j)*y_[j];
			}
		}

		for(unsigned int i(0);i<5;i++){
			std::cout<<y_[i]-extrapolate(c_[i])<<std::endl;
		}
	} else {
		std::cout<<"void Interpolation::compute_weights() : can't invert, support : "<<support_<<" (rcn="<<rcn<<")"<<std::endl;
	}
}

double Interpolation::extrapolate(Vector<double> const& x) const {
	double y(weights_.back());
	for(unsigned int i(0);i<N_;i++){
		y += weights_(i)*(this->*phi_)(sqrt((x-c_[i]).norm_squared())/support_);
	}
	for(unsigned int i(0);i<dim_;i++){ y += weights_(i+N_)*x(i); }
	return y;
}

void Interpolation::add_data(Vector<double> const& c, double const& y){
	bool add(true);
	for(unsigned int i(0);i<N_;i++){
		if(my::are_equal(c,c_[i])){ 
			add=false; i=N_; 
			std::cout<<"trying to add a value already entered"<<std::endl;
		}
	}
	if(add){
		c_.push_back(c);
		y_.push_back(y);
		N_++;
	}
}

void Interpolation::check_matrix_filling(){
	Rand<unsigned int> rnd(0,N_-1);
	unsigned int idx;
	unsigned int l;
	for(unsigned int i(0);i<10;i++){
		idx = rnd.get();
		l = 0;
		for(unsigned int j(0);j<N_;j++){
			if(sqrt((c_[idx]-c_[j]).norm_squared())<support_){
				l++;
			}
		}
		std::cout<<idx<<" "<<l<<std::endl;
	}
}

void Interpolation::add_data(unsigned int const& i, Vector<double> const& c, double const& y){
	c_[i]=c;
	y_[i]=y;
}

void Interpolation::set(unsigned int const& N, double const& linear_density){ 
	c_.clear();
	y_.clear();
	if(N){
		c_.resize(N);
		y_.resize(N);
		N_=N;
	} else { N_=0; }

	support_ = linear_density*pow(100*tgamma(dim_/2.+1)/pow(M_PI,dim_/2.),1./dim_);
	std::cout<<support_<<" "<<linear_density<<std::endl;
}

bool Interpolation::select_basis_function(unsigned int const& basis){
	switch(basis){
		case 1: { phi_ = &Interpolation::phi1; } break;
		case 2: { phi_ = &Interpolation::phi2; } break;
		case 3: { phi_ = &Interpolation::phi3; } break;
		case 4: { phi_ = &Interpolation::phi4; } break;
		case 5: { phi_ = &Interpolation::phi5; } break;
		case 6: { phi_ = &Interpolation::phi6; } break;
		case 7: { phi_ = &Interpolation::phi7; } break;
		case 8: { phi_ = &Interpolation::phi8; } break;
		case 9: { phi_ = &Interpolation::phi9; } break;
		default: {
					 std::cerr<<"Interpolation::Interpolation(unsigned int const& k) : unknown basis function"<<std::endl;
					 return false;
				 } break;
	}
	return true;
}
