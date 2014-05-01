#include "Binning.hpp"

/*constructors and destructor*/
/*{*/
Binning::Binning():
	B_(300),
	b_(5),
	bin_(new Vector<double>[b_]),
	log_(NULL)
{ set(); }

Binning::Binning(unsigned int B, unsigned int b):
	B_(B),
	b_(b),
	bin_(new Vector<double>[b_]),
	log_(NULL)
{ set(); }

void Binning::set(){
	x_ = 0.0;
	dx_ = 0.0;
	l_ = 0;
	DPL_ = 1;
	dpl_ = 0;
	log_iter_ = 0;
	Ml_.set(b_,0);
	m_bin_.set(b_,0);
	bin_[b_-1].set(2*B_,0);
	for(unsigned int i(b_-1);i>0;i--){
		bin_[i-1].set(2*bin_[i].size(),0);
	}
}

Binning::~Binning(){
	delete[] bin_;
	if(log_){ delete log_;}
}
/*}*/

/*public methods*/
/*{*/
void Binning::add_sample(double const& x){
	/*!add new entry to the bins*/
	if(DPL_ == ++dpl_){
		add_bin(0,2.0*x/DPL_,2.0*bin_[0](Ml_(0))/DPL_);
		dpl_ = 0;
		recompute_xdx_usefull_ = true;

		/*!update the bins if the bigger binning is big enough*/
		if(Ml_(b_-1)==2*B_){
			for(unsigned int i(0);i<b_-1;i++){
				for(unsigned int j(0);j<bin_[i+1].size();j++){
					bin_[i](j) = bin_[i+1](j);
				}
				for(unsigned int j(bin_[i+1].size());j<bin_[i].size();j++){
					bin_[i](j) = 0.0;/*I think it's useless */
				}
				Ml_(i) = Ml_(i+1);
				m_bin_(i) = m_bin_(i+1);
			}
			Ml_(b_-1)=0;
			m_bin_(b_-1)=0;
			for(unsigned int i(0);i<B_;i++){
				add_bin(b_-1,bin_[b_-2](2*i+1),bin_[b_-2](2*i));
			}
			for(unsigned int i(B_);i<2*B_;i++){
				bin_[b_-1](i) = 0;
			}
			l_++;
			DPL_*=2;
		}
	} else {
		bin_[0](Ml_(0)) += x; 
	}
}

bool Binning::is_converged(double const& tol){
	if(recompute_xdx_usefull_ && l_>0){
		/*!Compute the variance for each bin*/
		Vector<double> var_bin(b_,0.0);
		for(unsigned int l(0);l<b_;l++){
			for(unsigned int j(0);j<Ml_(l);j++){
				var_bin(l) += (bin_[l](j)-m_bin_(l))*(bin_[l](j)-m_bin_(l));
			}
			var_bin(l) = sqrt(var_bin(l) / (Ml_(l)*(Ml_(l)-1)));
		}
		/*!Do a linear regression to get the variance in the limit l->infty*/
		Vector<double> x(b_);
		double xb(0);
		double yb(0);
		for(unsigned int i(0);i<b_;i++){
			x(i) = 1.0/(l_+i);
			xb += x(i);
			yb += var_bin(i);
		}
		xb /= b_;
		yb /= b_;
		double num(0);
		double den(0);
		for(unsigned int i(0);i<b_;i++){
			num += (x(i)-xb)*(var_bin(i)-yb);
			den += (x(i)-xb)*(x(i)-xb);
		}

		x_=0.0;
		for(unsigned int i(0);i<Ml_(0);i++){ 
			x_ += bin_[0](i); 
		}
		x_/=Ml_(0);
		dx_ = yb-num/den*xb;
		recompute_xdx_usefull_ = false;

		if(log_){
			log_iter_++;
			(*log_)<<"$"+tostring(num/den)+"("+tostring(l_)<<")$"<<Write::endl;
			for(unsigned int i(0);i<b_;i++){(*log_)<<x(i)<<" "<<var_bin(i)<<Write::endl;}
			(*log_)<<Write::endl<<Write::endl;
		}

		if(std::abs(num/den)<tol){ return true; }
	} 
	return false;
}

void Binning::plot(){
	if(log_iter_>0){
		if(log_){
			Gnuplot gp("./",log_->get_filename(),"plot");
			gp.preplot("set xlabel '$\\ell^{-1}$'");
			gp.preplot("set ylabel '$\\Delta_\\ell$' rotate by 0 offset 2");
			gp.preplot("set xrange [0:1]");
			gp.preplot("set yrange [0:]");
			gp.add("for [IDX=0:"+tostring(log_iter_-1)+"] '"+log_->get_filename()+"' i IDX t columnheader(1)");
			//gp.add("'"+log_->get_filename()+"'");
			gp.save_file();
			gp.create_image(true);
		} else {
			std::cerr<<"Binning::plot() : no log file, can't create a plot"<<std::endl; 
		}
	} else {
		std::cerr<<"Binning::plot() : not enough sample to do a binning analysis"<<std::endl; 
	}
}
/*}*/

/*private methods*/
/*{*/
void Binning::add_bin(unsigned int l, double a, double b){
	recompute_xdx_usefull_ = true;
	bin_[l](Ml_(l)) = (a+b)/2.0;
	m_bin_(l) = (m_bin_(l)*Ml_(l)+bin_[l](Ml_(l)))/(Ml_(l)+1);
	Ml_(l)++;
	/*!Create the next (bigger) bin*/
	if(Ml_(l)%2==0 && l<b_-1){
		add_bin(l+1,bin_[l](Ml_(l)-1),bin_[l](Ml_(l)-2));
	}
}
/*}*/

