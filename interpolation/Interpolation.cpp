#include "Interpolation.hpp"

template<>
double Interpolation<double>::r(double const& a, double const& b) const { return std::abs(a-b)/support_; }

template<>
double Interpolation<Vector<double> >::r(Vector<double> const& a, Vector<double> const& b) const { return (a-b).norm_squared()/support_; }

template<>
double Interpolation<double>::operator()(double const& x) const {
	double y(0);
	if(weights_.size() == N_ +dim_+1){
		y = weights_.back();
		for(unsigned int i(0);i<N_;i++){
			y += weights_(i)*(this->*phi_)(r(x,c_[i]));
		}
		for(unsigned int i(0);i<dim_;i++){ y += weights_(i+N_)*x; }
	}
	if(N_ == weights_.size()){
		for(unsigned int i(0);i<N_;i++){
			y += weights_(i)*(this->*phi_)(r(x,c_[i]));
		}
	}
	return y;
}

template<>
void Interpolation<Vector<double> >::add_linear_term(Vector<double> const& x, double& y) const {
	for(unsigned int i(0);i<dim_;i++){ y += weights_(i+N_)*x(i); }
}

template<>
void Interpolation<double>::add_linear_term(double const& x, double& y) const {
	y += weights_(N_)*x;
}

template<>
void Interpolation<double>::add_linear_term(Matrix<double>& m){
	for(unsigned int i(0);i<N_;i++){ m(i,N_) = c_[i]; }
}

template<>
void Interpolation<Vector<double> >::add_linear_term(Matrix<double>& m){
	unsigned int k;
	for(unsigned int i(0);i<N_;i++){
		k=0;
		for(unsigned int j(N_);j<N_+dim_;j++){ m(i,j) = c_[i](k++); }
	}
}

template<>
void Interpolation<double>::add_linear_term(SparseMatrix<double>& m){
	for(unsigned int i(0);i<N_;i++){ m.push_back(i,N_,c_[i]); }
}

template<>
void Interpolation<Vector<double> >::add_linear_term(SparseMatrix<double>& m){
	unsigned int k;
	for(unsigned int i(0);i<N_;i++){
		k=0;
		for(unsigned int j(N_);j<N_+dim_;j++){ m.push_back(i,j,c_[i](k++)); }
	}
}
