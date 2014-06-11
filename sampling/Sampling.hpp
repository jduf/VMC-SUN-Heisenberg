#ifndef DEF_SAMPLING
#define DEF_SAMPLING

#include "Gnuplot.hpp"
#include "Vector.hpp"

template<typename Type>
class Binning;

template<typename Type>
class Binning{
	public:
		/*!Default constructor (needed for BinningSet)*/
		Binning(); 
		/*!Destructor*/
		~Binning();
		/*!Set the class*/
		void set(unsigned int const& B, unsigned int const& b);
		/*!Set a filename for the plot*/
		void plot(std::string const& filename){ log_ = new IOFiles(filename,true);}

		/*!Add sample to the bins*/
		void add_sample(Type const& x);
		/*!Compute the convergence*/
		void compute_convergence(double const& tol, Type& dx, bool& conv);
		/*!Set x_ to the mean value, dx_ to the variance, plot logl_ if any*/
		void complete_analysis(double const& tol, Type& x, Type& dx, bool& conv);
		/*!Compute the mean value*/
		double compute_x() const;

	private:
		unsigned int B_;	//!< minimum number of biggest bins needed to compute variance
		unsigned int b_;	//!< l_+b_ rank of the biggest bin (b !> 30)
		unsigned int l_;	//!< rank of the "smallest" bin 
		unsigned int DPL_;	//!< 2^l_ maximum number of element in each bin of rank l_
		unsigned int dpl_;	//!< current number of element in each bin of rank l_
		unsigned int logl_;	//!< number of Delta_l stored in log

		Vector<unsigned int> Ml_;//!<number bins of rank l : Ml = M0/2^l
		Vector<double> m_bin_;	//!< mean of the Binnings
		Vector<double>*bin_;	//!< Binnings

		bool recompute_dx_usefull_;	//!< true if dx should be recomputed
		bool addlog_;				//!< true if should be added to the log_

		IOFiles* log_;//!< log of dx_ when recompute_dx_usefull_

		/*!Recursive method that add samples in the different bins*/
		void add_bin(unsigned int l, double a, double b);
};

template<typename Type>
class Data{
	public:
		/*!Default constructor*/
		Data();
		/*!Destructor*/
		~Data();
		void set(unsigned int const& B, unsigned int const& b, bool const& conv);

		void add_sample(Data<Type> const& cs);
		void add_sample();

		void compute_convergence(double const& tol);
		void complete_analysis(double const& tol);
		void complete_analysis();

		/*maybe should return a const&*/
		Type const& get_x() const { return x_; } 
		Type const& get_dx() const { return dx_; } 
		bool const& get_conv() const { return conv_; } 
		unsigned int const& get_N() const { return N_; }

		void set_x(Type const& x){x_ = x;}
		void set_dx(Type const& dx){dx_ = dx;}
		void set_N(unsigned int const& N){N_ = N;}
		void set_conv(bool const& conv){conv_ = conv;}

		Data<Type>&  operator=(Type const& x){x_ =  x; return (*this);}
		Data<Type>& operator+=(Type const& x){x_ += x; return (*this);}
		Data<Type>& operator-=(Type const& x){x_ -= x; return (*this);}
		Data<Type>& operator/=(Type const& x){x_ /= x; return (*this);}
		Data<Type>& operator*=(Type const& x){x_ *= x; return (*this);}

		Binning<Type>* get_binning() const { return binning_; }

		void header_rst(std::string const& s, RST& rst) const;

	private:
		Data<Type>(Data<Type> const& d);

		Type x_;
		Type dx_;
		unsigned int N_;
		bool conv_; 
		Binning<Type>* binning_;
};

template<typename Type>
class DataSet{
	public:
		/*!Default constructor*/
		DataSet();
		/*!Destructor*/
		~DataSet();
		void set(unsigned int const& N);
		void set(unsigned int const& N, unsigned int const& B, unsigned int const& b, bool const& conv);

		void add_sample(DataSet<Type> const& ss);
		void add_sample();

		void compute_convergence(double const& tol);
		void complete_analysis(double const& tol);
		void complete_analysis();

		void header_rst(std::string const& s, RST& rst) const;

		unsigned int size() const { return size_;}
		
		Data<Type> const& operator[](unsigned int const& i) const 
		{assert(i<size_); return ds_[i];} 
		Data<Type>& operator[](unsigned int const& i) 
		{assert(i<size_); return ds_[i];} 

	private:
		Data<Type>* ds_;
		unsigned int size_;
};

/*Binning*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Binning<Type>::Binning():
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
void Binning<Type>::set(unsigned int const& B, unsigned int const& b){
	if(!bin_){
		B_ = B;
		b_ = b;
		bin_ = new Vector<double>[b];
	} else {
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

template<typename Type>
Binning<Type>::~Binning(){
	if(bin_){ delete[] bin_; }
	if(log_){ delete log_; }
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void Binning<Type>::add_sample(Type const& x){
	/*!add new entry to the bins*/
	if(DPL_ == ++dpl_){
		add_bin(0,2.0*x/DPL_,2.0*bin_[0](Ml_(0))/DPL_);
		recompute_dx_usefull_ = true;
		dpl_ = 0;
		/*!update the bins if the bigger Binning is big enough*/
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
		bin_[0](Ml_(0)) += x; 
	}
}

template<typename Type>
void Binning<Type>::compute_convergence(double const& tol, Type& dx, bool& conv){
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
		if(num/den>0.0){
			dx = yb - num/den * ( xb + 2.0*(l_+b_) ); 
		} else { dx = yb; }

		double criteria(num/den*compute_x());

		if(std::abs(criteria)<tol){ conv = true; }
		else{ conv = false; }

		if(log_ && addlog_){
			logl_++;
			addlog_ = false;
			(*log_)<<criteria<<IOFiles::endl;
			for(unsigned int i(0);i<b_;i++){(*log_)<<x(i)<<" "<<var_bin(i)<<IOFiles::endl;}
			(*log_)<<IOFiles::endl<<IOFiles::endl;
		}
	}
}

template<typename Type>
void Binning<Type>::complete_analysis(double const& tol, Type& x, Type& dx, bool& conv){
	recompute_dx_usefull_ = addlog_ = true;
	compute_convergence(tol,dx,conv);
	x = compute_x();
	if(l_>0 && log_){
		Gnuplot gp("./",log_->get_filename());
		gp.xrange(0,"");
		gp.yrange(0,"");
		gp+="set xlabel '$\\ell$'";
		gp+="set ylabel '$\\Delta_{\\ell}$' rotate by 0 offset 2";
		gp+="set key left bottom";
		gp+="plot for [IDX=0:"+tostring(logl_-1)+"] '"+log_->get_filename()+"' i IDX t columnheader(1)";
		if(log_){delete log_;}
		log_=NULL;
		gp.save_file();
		gp.create_image(true);
	} 
}

template<typename Type>
double Binning<Type>::compute_x() const {
	double x(0.0);
	for(unsigned int i(0);i<Ml_(0);i++){ x += bin_[0](i); }
	return x/Ml_(0); 
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void Binning<Type>::add_bin(unsigned int l, double a, double b){
	bin_[l](Ml_(l)) = (a+b)/2.0;
	m_bin_(l) = (m_bin_(l)*Ml_(l)+bin_[l](Ml_(l)))/(Ml_(l)+1);
	Ml_(l)++;
	/*!Create the next (bigger) bin*/
	if(Ml_(l)%2==0 && l<b_-1){
		add_bin(l+1,bin_[l](Ml_(l)-1),bin_[l](Ml_(l)-2));
	}
}
/*}*/
/*}*/

/*Data*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Data<Type>::Data():
	x_(0.0),
	dx_(0.0),
	N_(0),
	conv_(false),
	binning_(NULL)
{}

template<typename Type>
void Data<Type>::set(unsigned int const& B, unsigned int const& b, bool const& conv){
	binning_ = new Binning<Type>();
	binning_->set(B,b);
	conv_ = conv;
	N_ = 1;
}

template<typename Type>
Data<Type>::~Data(){
	if(binning_){ delete binning_; }
}
/*}*/

/*i/o methods*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, Data<Type> const& d){
	flux<<d.get_x()<<" "<<d.get_dx()<<" "<<d.get_N()<<" "<<d.get_conv(); 
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Data<Type>& d){
	Type x(0.0);
	Type dx(0.0);
	unsigned int N(0);
	bool conv(false);
	flux>>x>>dx>>N>>conv;
	d.set_x(x);
	d.set_dx(dx);
	d.set_N(N);
	d.set_conv(conv);
	return flux;
}

template<typename Type>
void Data<Type>::header_rst(std::string const& s, RST& rst) const {
	if(conv_){ rst.def(s,tostring(get_x())+" ("+tostring(get_dx())+")"); }
	else{ rst.def(s,"nc:"+tostring(get_x())+" ("+tostring(get_dx())+")"); }
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, Data<Type> const& d){
	if(w.is_binary()){ w<<d.get_x()<<d.get_dx()<<d.get_N()<<d.get_conv(); } 
	else { w.stream()<<d; }
	return w;
}

template<typename Type>
IOFiles& operator>>(IOFiles& r, Data<Type>& d){
	Type x(0.0);
	Type dx(0.0);
	unsigned int N(0);
	bool conv(false);
	r>>x>>dx>>N>>conv;
	d.set_x(x);
	d.set_dx(dx);
	d.set_N(N);
	d.set_conv(conv);
	return r;
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void Data<Type>::add_sample(Data<Type> const& cs){
	if(cs.conv_){
		x_ += cs.x_;
		dx_+= cs.dx_;
		N_ += cs.N_;
		conv_ = true;
	}
}

template<typename Type>
void Data<Type>::add_sample(){
	if(binning_){ binning_->add_sample(x_); }
}

template<typename Type>
void Data<Type>::compute_convergence(double const& tol) {
	if(binning_){ binning_->compute_convergence(tol,dx_,conv_);}
}

template<typename Type>
void Data<Type>::complete_analysis(double const& tol){
	if(binning_){ binning_->complete_analysis(tol,x_,dx_,conv_); }
}

template<typename Type>
void Data<Type>::complete_analysis(){
	x_ = x_/N_;
	dx_ = dx_/(N_*sqrt(N_));
}
/*}*/
/*}*/

/*DataSet*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
DataSet<Type>::DataSet():
	ds_(NULL),
	size_(0)
{}

template<typename Type>
void DataSet<Type>::set(unsigned int const& N){
	if(ds_){ delete[] ds_; }
	ds_ = new Data<Type>[N];
	size_ = N;
}

template<typename Type>
void DataSet<Type>::set(unsigned int const& N, unsigned int const& B, unsigned int const& b, bool const& conv){
	if(ds_){ delete[] ds_; }
	ds_ = new Data<Type>[N];
	size_ = N;
	for(unsigned int i(0);i<size_;i++){ ds_[i].set(B,b,conv); }
}

template<typename Type>
DataSet<Type>::~DataSet(){
	if(ds_){ delete[] ds_; }
}
/*}*/

/*i/o methods*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, DataSet<Type> const& ds){
	for(unsigned int i(0);i<ds.size();i++){ flux<<ds[i]<<std::endl; }
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, DataSet<Type> const& ds){
	std::cerr<<" std::istream& operator>>(std::istream& flux, DataSet<Type> const& v) not defined"<<std::endl;
	return flux;
}

template<typename Type>
void DataSet<Type>::header_rst(std::string const& s, RST& rst) const {
	rst.def(s,tostring(size_));
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, DataSet<Type> const& ds){
	if(w.is_binary()){
		w<<ds.size();
		for(unsigned int i(0);i<ds.size();i++){ w<<ds[i]; }
	} else {
		for(unsigned int i(0);i<ds.size();i++){ w.stream()<<ds[i]<<std::endl; }
	}
	return w;
}

template<typename Type>
IOFiles& operator>>(IOFiles& r, DataSet<Type>& ds){
	if(r.is_binary()){
		unsigned int size(0);
		r>>size;
		if(size != ds.size()){ ds.set(size); } 
		for(unsigned int i(0);i<ds.size();i++){ r>>ds[i]; }
	} else {
		for(unsigned int i(0);i<ds.size();i++){ r.stream()>>ds[i]; }
	}
	return r;
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void DataSet<Type>::add_sample(DataSet<Type> const& ds){
	assert(size_ == ds.size());
	for(unsigned int i(0);i<size_;i++){ ds_[i].add_sample(ds[i]); }
}

template<typename Type>
void DataSet<Type>::add_sample(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].add_sample(); }
}

template<typename Type>
void DataSet<Type>::compute_convergence(double const& tol){
	for(unsigned int i(0);i<size_;i++){ ds_[i].compute_convergence(tol); }
}

template<typename Type>
void DataSet<Type>::complete_analysis(double const& tol){
	for(unsigned int i(0);i<size_;i++){ ds_[i].complete_analysis(tol); }
}

template<typename Type>
void DataSet<Type>::complete_analysis(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].complete_analysis(); }
}
/*}*/
/*}*/
#endif
