#ifndef DEF_SAMPLING
#define DEF_SAMPLING

#include<cmath>
#include "Gnuplot.hpp"

template<typename Type>
class CorrelatedSamples;

template<typename Type>
class Data{
	public:
		/*!Default constructor*/
		Data();
		/*!Constructor*/
		//Data(Type const& x, Type const& dx, unsigned int const& N, bool const& conv=false);
		//void set(Type const& x, Type const& dx, unsigned int const& N, bool const& conv=false);

		void add_sample(CorrelatedSamples<Type> const& cs);

		/*maybe should return a const&*/
		Type get_x() const { return x_/N_;} 
		Type get_dx() const { return dx_/(N_*sqrt(N_));} /*to check formula and sqrt(N)=double*/
		bool get_conv() const { return conv_;} 

		void set_x(Type const& x){x_ = x;}
		void set_dx(Type const& dx){dx_ = dx;}
		void set_N(unsigned int const& N){N_ = N;}
		void set_conv(bool const& conv){conv_ = conv;}

	protected:
		Type x_;
		Type dx_;
		unsigned int N_;
		bool conv_; 
};

template<typename Type>
class CorrelatedSamples:public Data<Type>{
	public:
		/*!Default constructor (needed for CorrelatedSamplesSet)*/
		CorrelatedSamples(); 
		/*!Destructor*/
		~CorrelatedSamples();
		/*Set the whole class*/
		void set(unsigned int const& B, unsigned int const& b);

		void compute_convergence(double const& tol);

		/*Add sample to the bins*/
		void add_sample();

		void plot(std::string const& filename){ log_ = new Write(filename);}
		void complete_analysis(double const& tol);

		CorrelatedSamples<Type>&  operator=(Type const& x){this->x_ =  x; return (*this);}
		CorrelatedSamples<Type>& operator+=(Type const& x){this->x_ += x; return (*this);}
		CorrelatedSamples<Type>& operator-=(Type const& x){this->x_ -= x; return (*this);}
		CorrelatedSamples<Type>& operator/=(Type const& x){this->x_ /= x; return (*this);}
		CorrelatedSamples<Type>& operator*=(Type const& x){this->x_ *= x; return (*this);}

	private:
		unsigned int B_;	//!< minimum number of biggest bins needed to compute variance
		unsigned int b_;	//!< l_+b_ rank of the biggest bin (b !> 30)
		unsigned int l_;	//!< rank of the "smallest" bin 
		unsigned int DPL_;	//!< 2^l_ maximum number of element in each bin of rank l_
		unsigned int dpl_;	//!< current number of element in each bin of rank l_
		unsigned int logl_;	//!< number of Delta_l stored in log

		Vector<unsigned int> Ml_;//!<number bins of rank l : Ml = M0/2^l
		Vector<double> m_bin_;	//!< mean of the CorrelatedSampless
		Vector<double>*bin_;	//!< CorrelatedSampless

		bool recompute_dx_usefull_;	//!< true if dx should be recomputed
		bool addlog_;				//!< true if dx should be recomputed

		Write* log_;//!< log of dx_ when recompute_dx_usefull_

		/*!Recursive method that add samples in the different bins*/
		void add_bin(unsigned int l, double a, double b);
};

/*Data*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Data<Type>::Data():
	x_(0.0),
	dx_(0.0),
	N_(0),
	conv_(false)
{}

//template<typename Type>
//Data<Type>::Data(Type const& x, Type const& dx, unsigned int const& N, bool const& conv):
	//x_(x),
	//dx_(dx),
	//N_(N),
	//conv_(conv)
//{}

//template<typename Type>
//void Data<Type>::set(Type const& x, Type const& dx, unsigned int const& N, bool const& conv){
	//x_ = x;
	//dx_ = dx;
	//N_ = N;
	//conv_ = conv;
//}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void Data<Type>::add_sample(CorrelatedSamples<Type> const& cs){
	if(this->conv_ == cs.conv_){
		this->x_ += cs.x_;
		this->dx_+= cs.dx_;
		this->N_ += cs.N_;
	}
}
/*}*/
/*}*/

/*CorrelationSamples*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
CorrelatedSamples<Type>::CorrelatedSamples():
	B_(0),
	b_(0),
	l_(0),
	DPL_(1),
	dpl_(0),
	logl_(0),
	bin_(NULL),
	recompute_dx_usefull_(false),
	addlog_(false),
	log_(NULL)
{}

template<typename Type>
CorrelatedSamples<Type>::~CorrelatedSamples(){
	if(bin_){delete[] bin_;}
	if(log_){delete log_;}
}

template<typename Type>
void CorrelatedSamples<Type>::set(unsigned int const& B, unsigned int const& b){
	if(!bin_){
		B_ = B;
		b_ = b;
		bin_ = new Vector<double>[b];
	} else {
		this->x_ = 0.0;
		this->dx_ = 0.0;
		this->N_ = 0;
		this->conv_ = false;
		l_ = 0;
		DPL_ = 1;
		dpl_ = 0;
		logl_ = 0;
		recompute_dx_usefull_ = false;
		addlog_ = false;
	}
	Ml_.set(b_,0);
	m_bin_.set(b_,0);
	bin_[b_-1].set(2*B_,0);
	for(unsigned int i(b_-1);i>0;i--){
		bin_[i-1].set(2*bin_[i].size(),0);
	}
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void CorrelatedSamples<Type>::add_sample(){
	/*!add new entry to the bins*/
	if(DPL_ == ++dpl_){
		add_bin(0,2.0*this->x_/DPL_,2.0*bin_[0](Ml_(0))/DPL_);
		dpl_ = 0;
		/*!update the bins if the bigger CorrelatedSamples is big enough*/
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
			addlog_ = true;
		}
	} else {
		bin_[0](Ml_(0)) += this->x_; 
	}
}

template<typename Type>
void CorrelatedSamples<Type>::compute_convergence(double const& tol){
	if(recompute_dx_usefull_ && Ml_(b_-1)>B_){
		recompute_dx_usefull_ = false;
		/*!Compute the variance for each bin*/
		Vector<double> var_bin(b_,0.0);
		for(unsigned int l(0);l<b_;l++){
			for(unsigned int j(0);j<Ml_(l);j++){
				var_bin(l) += (bin_[l](j)-m_bin_(l))*(bin_[l](j)-m_bin_(l));
			}
			var_bin(l) = sqrt(var_bin(l) / (Ml_(l)*(Ml_(l)-1)));
		}

		/*!Do a linear regression*/
		Vector<double> x(b_);
		double xb(0);
		double yb(0);
		double ymax(var_bin.max());;
		for(unsigned int i(0);i<b_;i++){
			x(i) = l_+i;
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
		double tmp(yb - num/den * xb);

		if(num/den>0.0){ this->dx_ = tmp + num/den*2.0*(l_+b_); } 
		else { this->dx_ = yb; }
		if(std::abs(num/(den*ymax))<tol){ this->conv_ = true; }

		if(log_ && addlog_){
			logl_++;
			addlog_ = false;
			(*log_)<<num/(den*ymax)<<Write::endl;
			for(unsigned int i(0);i<b_;i++){(*log_)<<x(i)<<" "<<var_bin(i)/ymax<<Write::endl;}
			(*log_)<<Write::endl<<Write::endl;
		}
	}
}

template<typename Type>
void CorrelatedSamples<Type>::complete_analysis(double const& tol){
	compute_convergence(tol);
	this->x_ = 0.0;
	for(unsigned int i(0);i<Ml_(0);i++){ this->x_ += bin_[0](i); }
	this->x_ /= Ml_(0); 
	this->N_ = 1;
	if(l_>0 && log_){
		Gnuplot gp("./",log_->get_filename());
		gp.xrange(0,"");
		gp.yrange(0,1.1);
		gp+="set xlabel '$\\ell^{-1}$'";
		gp+="set ylabel '$\\Delta_\\ell$' rotate by 0 offset 2";
		gp+="set key left bottom";
		gp+="plot for [IDX=0:"+tostring(logl_-1)+"] '"+log_->get_filename()+"' i IDX t columnheader(1)";
		gp.save_file();
		gp.create_image(true);
	} 
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void CorrelatedSamples<Type>::add_bin(unsigned int l, double a, double b){
	bin_[l](Ml_(l)) = (a+b)/2.0;
	m_bin_(l) = (m_bin_(l)*Ml_(l)+bin_[l](Ml_(l)))/(Ml_(l)+1);
	Ml_(l)++;
	/*!Create the next (bigger) bin*/
	if(Ml_(l)%2==0 && l<b_-1){
		add_bin(l+1,bin_[l](Ml_(l)-1),bin_[l](Ml_(l)-2));
		recompute_dx_usefull_ = true;
	}
}
/*}*/
/*}*/
#endif
