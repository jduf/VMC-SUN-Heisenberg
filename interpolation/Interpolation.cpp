#include "Interpolation.hpp"

Interpolation::Interpolation(unsigned int const& dim):
	basis_(0),
	dim_(dim),
	N_(0),
	support_(0),
	phi_(NULL)
{}

void Interpolation::set_data(){
	c_.clear();
	y_.clear();
	N_ = 0;
}

bool Interpolation::compute_weights(double const& dx, unsigned int const& n){
	if(basis_>6){
		SparseMatrix<double> m;
		support_ = dx*n*pow(50*tgamma(dim_/2.+1)/pow(M_PI,dim_/2.)/N_,1./dim_);
		double radius;
		unsigned int k;
		for(unsigned int i=0;i<N_;i++){
			m.push_back(i,i,1);
			for(unsigned int j(i+1);j<N_;j++){
				radius = sqrt((c_[i]-c_[j]).norm_squared())/support_;
				if(radius<1){ m.push_back(i,j,(this->*phi_)(radius)); }
			}
			k=0;
			for(unsigned int j(N_);j<N_+dim_;j++){ m.push_back(i,j,c_[i](k++)); }
			m.push_back(i,N_+dim_,1.0);
		}

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
			if(row!=col){r0(col) -= m[idx]*weights_(row);}
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

		double err(0);
		unsigned int idx;
		unsigned int n_err(100);
		Rand<unsigned int> rnd_idx(0,N_-1);
		for(unsigned int i(0);i<n_err;i++){
			idx = rnd_idx.get();
			err += std::abs(y_[idx]-extrapolate(c_[idx]));
		}
		err /= n_err;
		if(err>1e-14){
			std::cerr<<"void Interpolation::compute_weights(double const& dx, unsigned int const& n) : warning error : "<<err<<std::endl;
		}
		return true;
	} else {
		std::cout<<"void Interpolation::compute_weights(double const& dx, unsigned int const& n) : not a CBRF"<<std::endl;
		return false;
	}
}

bool Interpolation::compute_weights(){
	if(basis_<=6){
		Vector<int> ipiv;
		Matrix<double> m(N_+dim_+1,N_+dim_+1);
		double rcn;
		unsigned int k(0);
		support_ = 1.0;
#pragma omp parallel for private(k)
		for(unsigned int i=0;i<N_;i++){
			for(unsigned int j(i+1);j<N_;j++){
				m(i,j) = (this->*phi_)(sqrt((c_[i]-c_[j]).norm_squared()));
			}
			k=0;
			for(unsigned int j(N_);j<N_+dim_;j++){ m(i,j) = c_[i](k++); }
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
				err += std::abs(y_[idx]-extrapolate(c_[idx]));
			}
			err /= n_err;
			if(err>1e-14){
				std::cerr<<"void Interpolation::compute_weights() : warning error : "<<err<<std::endl;
			}
			return true;
		} else { std::cout<<"void Interpolation::compute_weights() : can't invert (rcn="<<rcn<<")"<<std::endl; }
	} else { std::cout<<"void Interpolation::compute_weights() : not a BRF, needs a support"<<std::endl; }
	return false;
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

bool Interpolation::select_basis_function(unsigned int const& basis){
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
					 std::cerr<<"Interpolation::Interpolation(unsigned int const& k) : unknown basis function"<<std::endl;
					 return false;
				 } break;
	}
	return true;
}
