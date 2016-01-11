#ifndef DEF_INTERPOLATION
#define DEF_INTERPOLATION

#include "Lapack.hpp"
#include "omp.h"
#include "Rand.hpp"

/*!If two data <c,y> are identical the matrix to invert won't be invertible*/
template<typename Type>
class Interpolation{
	public:
		Interpolation(unsigned int const& dim);

		bool select_basis_function(unsigned int const& basis);

		void add_data(Type const& ci, double const& yi);
		void add_data(unsigned int const& i, Type const& c, double const& y);

		void set_data();
		bool compute_weights();
		bool compute_weights(double& dx, unsigned int const& n);

		double operator()(Type const& x) const;

		unsigned int const& get_N() const { return N_; }

	private:
		std::vector<Type> c_;
		std::vector<double> y_;
		Vector<double> weights_;
		unsigned int basis_;
		unsigned int dim_;
		unsigned int N_;
		double support_;
		double (Interpolation::*phi_)(double const&) const;

		template<typename MType>
			class SparseMatrix{
				public:
					SparseMatrix() = default;

					void get_idx(unsigned int const& idx, unsigned int& r, unsigned int& c) const { r=r_[idx]; c=c_[idx]; }
					void set_idx(unsigned int const& idx, unsigned int const& r, unsigned int const& c, MType const& v){ r_[idx]=r; c_[idx]=c; v_[idx]=v; }

					void push_back(unsigned int const& r, unsigned int const& c, MType v){ r_.push_back(r); c_.push_back(c); v_.push_back(v); }

					MType const& operator[](unsigned int const& idx) const { return v_[idx]; }
					MType& operator[](unsigned int const& idx){ return v_[idx]; }

					unsigned int size() const { return v_.size(); }

				protected:
					std::vector<MType> v_;
					std::vector<unsigned int> r_;
					std::vector<unsigned int> c_;
			};

		double r(Type const& a, Type const& b) const;
		void add_linear_term(Type const& x, double& y) const;
		void add_linear_term(SparseMatrix<double>& m);
		void add_linear_term(Matrix<double>& m);

		double phi1(double const& r) const { return r; }
		double phi2(double const& r) const { return r*(r<1?log(pow(r,r)):r*log(r)); }
		double phi3(double const& r) const { return r*r*r; }
		double phi4(double const& r) const { return r*r*r*(r<1?log(pow(r,r)):r*log(r)); }
		double phi5(double const& r) const { return r*r*r*r*r; }
		double phi6(double const& r) const { return r*r*r*r*r*(r<1?log(pow(r,r)):r*log(r)); }
		double phi7(double const& r) const { return (r<1.0?(1-r)*(1-r):0.0); }
		double phi8(double const& r) const { return (r<1.0?(1-r)*(1-r)*(1-r)*(1-r)*(4*r+1):0.0); }
		double phi9(double const& r) const { return (r<1.0?(1-r)*(1-r)*(1-r)*(1-r)*(1-r)*(1-r)*(32*r*r*r+25*r*r+8*r+1):0.0); }
};

template<typename Type>
Interpolation<Type>::Interpolation(unsigned int const& dim):
	basis_(0),
	dim_(dim),
	N_(0),
	support_(0),
	phi_(NULL)
{}

template<typename Type>
void Interpolation<Type>::set_data(){
	c_.clear();
	y_.clear();
	N_ = 0;
}

template<typename Type>
void Interpolation<Type>::add_data(Type const& c, double const& y){
	bool add(true);
	for(unsigned int i(0);i<N_;i++){
		if(my::are_equal(c,c_[i])){
			std::cout<<"trying to add a value already entered "<<c<<" <-> "<<c_[i]<<std::endl;
			add=false;
			i=N_;
		}
	}
	if(add){
		c_.push_back(c);
		y_.push_back(y);
		N_++;
	}
}

template<typename Type>
void Interpolation<Type>::add_data(unsigned int const& i, Type const& c, double const& y){
	c_[i]=c;
	y_[i]=y;
}

template<typename Type>
bool Interpolation<Type>::compute_weights(double& dx, unsigned int const& n){
	if(basis_>6){
		SparseMatrix<double> m;
		support_ = dx*n*pow(50*tgamma(dim_/2.+1)/pow(M_PI,dim_/2.)/N_,1./dim_);
		double radius;
		for(unsigned int i=0;i<N_;i++){
			m.push_back(i,i,1);
			for(unsigned int j(i+1);j<N_;j++){
				radius = r(c_[i],c_[j]);
				if(radius<1){ m.push_back(i,j,(this->*phi_)(radius)); }
			}
			m.push_back(i,N_+dim_,1.0);
		}
		add_linear_term(m);

		Vector<double> r0(N_+dim_+1);
		for(unsigned int i(0);i<N_;i++){ r0(i) = y_[i]; }
		for(unsigned int i(N_);i<N_+dim_+1;i++){ r0(i) = 0; }
		Vector<double> p0;
		Vector<double> Ap;
		weights_.set(N_+dim_+1);
		unsigned int maxiter(1e5);
		unsigned int row;
		unsigned int col;
		double alpha;
		double beta;
		double r0_ns;
		double r1_ns;
		Rand<double> rnd(-1.0,1.0);

		for(unsigned int i(0);i<N_+dim_+1;i++){ weights_(i) = rnd.get(); }
		for(unsigned int idx(0);idx<m.size();idx++){
			m.get_idx(idx,row,col);
			r0(row) -= m[idx]*weights_(col);
			if(row!=col){ r0(col) -= m[idx]*weights_(row); }
		}
		p0 = r0;
		r0_ns = r0.norm_squared();

		for(unsigned int iter(0);iter<maxiter;iter++){
			Ap.set(N_+dim_+1,0);
			for(unsigned int idx(0);idx<m.size();idx++){
				m.get_idx(idx,row,col);
				Ap(row) += m[idx]*p0(col);
				if(row!=col){ Ap(col) += m[idx]*p0(row); }
			}

			alpha = 0;
			for(unsigned int i(0);i<N_+dim_+1;i++){ alpha += p0(i)*Ap(i); }
			alpha = r0_ns/alpha;
			weights_ += p0*alpha;
			r0 -= Ap*alpha;

			if( r0_ns<1e-40){ iter=maxiter; }
			else {
				r1_ns = r0.norm_squared();
				beta = r1_ns/r0_ns;
				r0_ns = r1_ns;

				p0 *= beta;
				p0 += r0;
			}
		}

		dx = 0;
		unsigned int idx;
		unsigned int n_err(100);
		Rand<unsigned int> rnd_idx(0,N_-1);
		for(unsigned int i(0);i<n_err;i++){
			idx = rnd_idx.get();
			dx += std::abs(y_[idx]-(*this)(c_[idx]));
		}
		dx /= n_err;
		if(dx>1e-14){
			std::cerr<<__PRETTY_FUNCTION__<<" : warning error : "<<dx<<std::endl;
		}
		return true;
	} else {
		std::cout<<"void Interpolation::compute_weights(double const& dx, unsigned int const& n) : not a CBRF"<<std::endl;
		return false;
	}
}

template<typename Type>
bool Interpolation<Type>::compute_weights(){
	support_ = 1.0;
	if(basis_<=6){
		Matrix<double> m(N_+dim_+1,N_+dim_+1);
#pragma omp parallel for
		for(unsigned int i=0;i<N_;i++){
			for(unsigned int j(i+1);j<N_;j++){
				m(i,j) = (this->*phi_)(r(c_[i],c_[j]));
			}
			m(i,N_+dim_) = 1.0;
		}
		add_linear_term(m);
		for(unsigned int i(N_);i<N_+dim_+1;i++){
			for(unsigned int j(i);j<N_+dim_+1;j++){
				m(i,j) = 0.0;
			}
		}

		Lapack<double> inv_m(m,false,'S');
		double rcn;
		Vector<int> ipiv(inv_m.is_singular(rcn));

		if(ipiv.ptr()){
			inv_m.inv(ipiv);
			weights_.set(N_+dim_+1,0);
#pragma omp parallel for
			for(unsigned int i=0;i<N_+dim_+1;i++){
				for(unsigned int j(0);j<N_;j++){
					weights_(i) += m(i,j)*y_[j];
				}
			}

			double err(0);
			unsigned int idx;
			unsigned int n_err(100);
			Rand<unsigned int> rnd_idx(0,N_-1);
			for(unsigned int i(0);i<n_err;i++){
				idx = rnd_idx.get();
				err += std::abs(y_[idx]-(*this)(c_[idx]));
			}
			err /= n_err;
			if(err>1e-14){
				std::cerr<<__PRETTY_FUNCTION__<<" : warning error : "<<err<<std::endl;
			}
			return true;
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : can't invert (rcn="<<rcn<<")"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<": not a BRF, needs a support"<<std::endl; }
	return false;
}

template<typename Type>
bool Interpolation<Type>::select_basis_function(unsigned int const& basis){
	basis_ = basis;
	switch(basis_){
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
					 std::cerr<<__PRETTY_FUNCTION__<<" : unknown basis function"<<std::endl;
					 return false;
				 } break;
	}
	return true;
}

template<typename Type>
double Interpolation<Type>::operator()(Type const& x) const {
	double y(0);
	y = weights_.back();
	for(unsigned int i(0);i<N_;i++){
		y += weights_(i)*(this->*phi_)(r(x,c_[i]));
	}
	if(weights_.size() == N_ +dim_+1){ add_linear_term(x,y); }
	return y;
}
#endif
