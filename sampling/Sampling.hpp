#ifndef DEF_SAMPLING
#define DEF_SAMPLING

#include "Vector.hpp"

/*{Description*/
/*!The use of this class with int or unsigned int is problemeatic because
 * during the binning procedure, an operation (a+b)/2 will have rounding
 * issues
 */
/*}*/
template<typename Type>
class Binning{
	public:
		/*!Constructor*/
		Binning(unsigned int const& B, unsigned int const& b); 
		/*!Copy constructor*/
		Binning(Binning const& b); 
		/*!Move constructor*/
		Binning(Binning&& b); 
		/*!Constructor that reads from file*/
		Binning(IOFiles& r);
		/*!Destructor*/
		~Binning();
		/*{Forbidden*/
		Binning& operator=(Binning) = delete;
		/*}*/

		/*!Set the class*/
		void set();

		/*!Add sample to the bins*/
		void add_sample(Type const& x);
		/*!Compute the convergence*/
		void compute_convergence(double const& tol, Type& dx, bool& conv);
		/*!Set x_ to the mean value, dx_ to the variance*/
		void complete_analysis(double const& tol, Type& x, Type& dx, bool& conv);
		/*!Compute the mean value*/
		Type const& get_x() const { return m_bin_(0); }
		/*!Merge this with b*/
		void merge(Binning const& b);

		void write(IOFiles& w) const;

	private:
		unsigned int const B_;//!< minimum number of biggest bins needed to compute variance
		unsigned int const b_;//!< l_+b_ rank of the biggest bin (b !> 30)
		unsigned int l_;	//!< rank of the "smallest" bin 
		unsigned int DPL_;	//!< 2^l_ maximum number of element in each bin of rank l_
		unsigned int dpl_;	//!< current number of element in each bin of rank l_

		Vector<unsigned int> Ml_;//!<number bins of rank l : Ml = M0/2^l
		Vector<Type> m_bin_;	//!< mean of the Binnings
		Vector<Type>*  bin_;	//!< Binnings

		bool recompute_dx_usefull_;	//!< true if dx should be recomputed

		/*!Recursive method that add samples in the different bins*/
		void add_bin(unsigned int l, Type const& a, Type const& b);
		void do_merge(Vector<Type> const& bin, unsigned int const& dpl, unsigned int const& DPL, unsigned int const& Ml);
};

template<typename Type>
class Data{
	public:
		/*!Default constructor*/
		Data();
		/*!Copy constructor*/
		Data<Type>(Data<Type> const& d);
		/*!Move constructor*/
		Data<Type>(Data<Type>&& d);
		/*!Constructor that reads from file*/
		Data<Type>(IOFiles& r);
		/*!Destructor*/
		~Data();
		void set(Type const& x, Type const& dx, unsigned int const& N, bool const& conv);
		void set(unsigned int const& B, unsigned int const& b, bool const& conv);
		void set();
		Data<Type>& operator=(Data<Type> d);

		void add_sample(Data<Type> const& d);
		void add_sample();
		void merge(Data const& d);

		void compute_convergence(double const& tol);
		void complete_analysis(double const& tol);
		void complete_analysis();
		void delete_binning();

		Type const& get_x() const { return binning_?binning_->get_x():x_; } 
		Type const& get_dx() const { return dx_; } 
		bool const& get_conv() const { return conv_; } 
		unsigned int const& get_N() const { return N_; }

		void set_x(Type const& x){x_ = x;}

		void add(Type const& x){x_ += x;}
		void substract(Type const& x){ x_ -= x; }
		void multiply(Type const& x){ 
			x_ *= x;
			if(!binning_){ dx_ *= sqrt(x); }
		}
		void divide(Type const& x){ 
			x_ /= x; 
			if(!binning_){ dx_ /= sqrt(x); }
		}

		Binning<Type>* get_binning() const { return binning_; }
		void set_binning(Binning<Type>* b) { binning_ = b; }

		void write(IOFiles& w, std::string const& name="") const;

		void header_rst(std::string const& s, RST& rst) const;

	private:
		Type x_;
		Type dx_;
		unsigned int N_;
		bool conv_; 
		Binning<Type>* binning_;

		void swap_to_assign(Data<Type>& d1, Data<Type>& d2);
};

template<typename Type>
class DataSet{
	public:
		/*!Default constructor*/
		DataSet();
		/*!Copy constructor*/
		DataSet<Type>(DataSet<Type> const& d);
		/*!Copy constructor*/
		DataSet<Type>(DataSet<Type>&& d);
		/*!Constructor that reads from file*/
		DataSet<Type>(IOFiles& r);
		/*!Destructor*/
		~DataSet();
		/*{Forbidden*/
		DataSet<Type>& operator=(DataSet<Type>) = delete;
		/*}*/

		void set();
		void set(unsigned int const& N);
		void set(unsigned int const& N, unsigned int const& B, unsigned int const& b, bool const& conv);

		void add_sample(DataSet<Type> const& ds);
		void add_sample();
		void merge(DataSet<Type> const& ds);

		void compute_convergence(double const& tol);
		void complete_analysis(double const& tol);
		void complete_analysis();
		void delete_binning();

		void header_rst(std::string const& s, RST& rst) const;

		unsigned int size() const { return size_;}

		Data<Type> const& operator[](unsigned int const& i) const 
		{assert(i<size_); return ds_[i];} 
		Data<Type>& operator[](unsigned int const& i) 
		{assert(i<size_); return ds_[i];} 

	private:
		unsigned int size_;
		Data<Type>* ds_;
};

/*Binning*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Binning<Type>::Binning(unsigned int const& B, unsigned int const& b):
	B_(B),
	b_(b),
	bin_(new Vector<Type>[b])
{set();}

template<typename Type>
Binning<Type>::Binning(Binning const& b):
	B_(b.B_),
	b_(b.b_),
	l_(b.l_),
	DPL_(b.DPL_),
	dpl_(b.dpl_),
	Ml_(b.Ml_),
	m_bin_(b.m_bin_),
	bin_(b.bin_?new Vector<Type>[b_]:NULL),
	recompute_dx_usefull_(b.recompute_dx_usefull_)
{
	for(unsigned int i(0);i<b_;i++){ bin_[i] = b.bin_[i]; }
	std::cout<<"binning copy"<<std::endl;
}

template<typename Type>
Binning<Type>::Binning(Binning&& b):
	B_(b.B_),
	b_(b.b_),
	l_(b.l_),
	DPL_(b.DPL_),
	dpl_(b.dpl_),
	Ml_(b.Ml_),
	m_bin_(std::move(b.m_bin_)),
	bin_(b.bin_),
	recompute_dx_usefull_(b.recompute_dx_usefull_)
{
	b.bin_ = NULL;
}

template<typename Type>
Binning<Type>::Binning(IOFiles& r):
	B_(r.read<unsigned int>()),
	b_(r.read<unsigned int>()),
	l_(r.read<unsigned int>()),
	DPL_(r.read<unsigned int>()),
	dpl_(r.read<unsigned int>()),
	Ml_(r),
	m_bin_(r),
	bin_(b_?new Vector<Type>[b_]:NULL),
	recompute_dx_usefull_(true)
{
	r>>bin_[0];
	for(unsigned int l(0);l<b_-1;l++){
		bin_[l+1].set(bin_[l].size()/2);
		for(unsigned int i(0);i<Ml_(l);i+=2){
			bin_[l+1](i/2) = (bin_[l](i)+bin_[l](i+1))/2.0;
		}
	}
}

template<typename Type>
void Binning<Type>::set(){
	l_ = 0;
	DPL_ = 1;
	dpl_ = 0;
	recompute_dx_usefull_ = false;
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
}
/*}*/

/*public methods*/
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
			for(unsigned int i(B_);i<2*B_;i++){ bin_[b_-1](i) = 0; }
			l_++;
			DPL_*=2;
		}
	} else { bin_[0](Ml_(0)) += x; }
}

template<typename Type>
void Binning<Type>::merge(Binning const& other){
	if(B_ == other.B_ && b_ == other.b_ ){
		if(l_>= other.l_){
			do_merge(other.bin_[0],other.dpl_,other.DPL_,other.Ml_(0));
		} else {
			Vector<Type> tmp_bin(bin_[0]);
			unsigned int tmp_Ml(Ml_(0));
			unsigned int tmp_dpl(dpl_);
			unsigned int tmp_DPL(DPL_);
			for(unsigned int i(0);i<b_;i++){ bin_[i] = other.bin_[i]; }
			Ml_ = other.Ml_;
			DPL_ = other.DPL_;
			l_ = other.l_;
			recompute_dx_usefull_ = other.recompute_dx_usefull_;
			do_merge(tmp_bin,tmp_dpl,tmp_DPL,tmp_Ml);
		}
	} else {
		std::cerr<<"void Binning<Type>::merge(Binning const& b) : B_ != b.B_ || b_ != b.b_ "<<std::endl;
	}
}

template<typename Type>
void Binning<Type>::compute_convergence(double const& tol, Type& dx, bool& conv){
	if(recompute_dx_usefull_ && Ml_(b_-1)>=B_){
		recompute_dx_usefull_ = false;
		/*!Compute the variance for each bin*/
		Vector<Type> var_bin(b_,0.0);
		for(unsigned int l(0);l<b_;l++){
			for(unsigned int j(0);j<Ml_(l);j++){
				var_bin(l) += my::norm_squared(bin_[l](j)-m_bin_(l));
			}
			var_bin(l) = sqrt(var_bin(l) / (Ml_(l)*(Ml_(l)-1)));
		}

		/*!Do a linear regression and if the slope is almost flat, the system
		 * is believed to be converged*/
		Vector<Type> x(b_);
		Type xb(0);
		Type yb(0);
		Type num(0);
		Type den(0);
		for(unsigned int i(0);i<b_;i++){
			x(i) = l_+i;
			xb += x(i);
			yb += var_bin(i);
		}
		xb /= b_;
		yb /= b_;
		for(unsigned int i(0);i<b_;i++){
			num += (x(i)-xb)*(var_bin(i)-yb);
			den += (x(i)-xb)*(x(i)-xb);
		}
		dx = yb;
		if(num/den>0.0){ dx += num/den * xb; } 
		Type criteria(num/den*m_bin_(0));
		if(std::abs(criteria)<tol){ conv = true; }
		else{ conv = false; }
	}
}

template<typename Type>
void Binning<Type>::complete_analysis(double const& tol, Type& x, Type& dx, bool& conv){
	recompute_dx_usefull_ = true;
	compute_convergence(tol,dx,conv);
	/*{old way to set x*/
	/*! x = m_bin_(0); */
	/*}*/
	x = (m_bin_(0)*Ml_(0)*DPL_+bin_[0](Ml_(0)))/(Ml_(0)*DPL_+dpl_);
	/*{Description*/
	/*!
	  std::cout<<"given x"<<x<<std::endl;
	  std::cout<<l_<<"leftover "<<bin_[0](Ml_(0))<<" with "<<dpl_<<" therefore "<<(x*Ml_(0)*DPL_+bin_[0](Ml_(0)))/(Ml_(0)*DPL_+dpl_)<<std::endl;
	  Type tmp(0);
	  for(unsigned int i(0);i<Ml_(0);i++){
	  tmp += bin_[0](i);
	  }
	  std::cout<<tmp<<" "<<tmp/Ml_(0)<<" "<<Ml_(0)<<std::endl;
	  */
	/*}*/
}

template<typename Type>
void Binning<Type>::write(IOFiles& w) const {
	w<<B_<<b_<<l_<<DPL_<<dpl_<<Ml_<<m_bin_<<bin_[0];
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void Binning<Type>::add_bin(unsigned int l, Type const& a, Type const& b){
	bin_[l](Ml_(l)) = (a+b)/2.0;
	m_bin_(l) = (m_bin_(l)*Ml_(l)+bin_[l](Ml_(l)))/(Ml_(l)+1);
	Ml_(l)++;
	/*!Create the next (bigger) bin*/
	if(Ml_(l)%2==0 && l<b_-1){
		add_bin(l+1,bin_[l](Ml_(l)-1),bin_[l](Ml_(l)-2));
	}
}

template<typename Type>
void Binning<Type>::do_merge(Vector<Type> const& bin, unsigned int const& dpl, unsigned int const& DPL, unsigned int const& Ml){
	Type old_last_bin(bin_[0](Ml_(0)));
	unsigned int old_dpl(dpl_);
	/*{Description*/
	/*
	   std::cout<<"sould be an integer "<<1.0*DPL_/(1.0*DPL)<<std::endl;
	   std::cout<<"CURRENT"<<std::endl;
	   std::cout<<"x    = "<<old_last_bin<<std::endl;
	   std::cout<<"dpl_ = "<<old_dpl<<std::endl;
	   std::cout<<"DPL_ = "<<DPL_<<std::endl;
	   std::cout<<"last = "<<Ml_(0)<<std::endl;
	   std::cout<<"OTHER"<<std::endl;
	   std::cout<<"x    = "<<bin(Ml)<<std::endl;
	   std::cout<<"dpl_ = "<<dpl<<std::endl;
	   std::cout<<"DPL_ = "<<DPL<<std::endl;
	   std::cout<<"last = "<<Ml<<std::endl;
	   */
	/*}*/
	bin_[0](Ml_(0)) = 0;
	dpl_ = 0;
	for(unsigned int i(0);i<Ml;i++){
		dpl_ += DPL-1;
		add_sample(DPL*bin(i));
	}
	Type tmp((bin(Ml)+old_last_bin)/(dpl+old_dpl));
	for(unsigned int i(0);i<dpl+old_dpl;i++){ add_sample(tmp); }
	/*{Description*/
	/*
	   std::cout<<"FINAL"<<std::endl;
	   std::cout<<"x    = "<<bin_[0](Ml_(0))<<std::endl;
	   std::cout<<"dpl_ = "<<dpl_<<std::endl;
	   std::cout<<"DPL_ = "<<DPL_<<std::endl;
	   std::cout<<"last = "<<Ml_(0)<<std::endl;
	   */
	/*}*/
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
Data<Type>::Data(Data<Type> const& d):
	x_(d.x_),
	dx_(d.dx_),
	N_(d.N_),
	conv_(d.conv_),
	binning_(d.binning_?new Binning<Type>(*d.binning_):NULL)
{
	std::cout<<"data copy"<<std::endl;
}

template<typename Type>
Data<Type>::Data(Data<Type>&& d):
	x_(d.x_),
	dx_(d.dx_),
	N_(d.N_),
	conv_(d.conv_),
	binning_(d.binning_)
{
	d.binning_ = NULL;
}

template<typename Type>
Data<Type>::Data(IOFiles& r):
	x_(r.read<Type>()),
	dx_(r.read<Type>()),
	N_(r.read<unsigned int>()),
	conv_(r.read<bool>()),
	binning_(r.read<bool>()?new Binning<Type>(r):NULL)
{}

template<typename Type>
void Data<Type>::set(Type const& x, Type const& dx, unsigned int const& N, bool const& conv){
	x_ = x;
	dx_ = dx;
	N_ = N;
	conv_ = conv;
}

template<typename Type>
void Data<Type>::set(unsigned int const& B, unsigned int const& b, bool const& conv){
	if(binning_){ delete binning_; }
	binning_ = new Binning<Type>(B,b);
	conv_ = conv;
	x_ = 0.0;
	dx_ = 0.0;
	N_ = 1;
}

template<typename Type>
void Data<Type>::set(){
	binning_->set();
	x_ = 0.0;
	dx_ = 0.0;
	N_ = 1;
}

template<typename Type>
Data<Type>::~Data(){
	if(binning_){ delete binning_; }
}

template<typename Type>
void Data<Type>::swap_to_assign(Data<Type>& d1, Data<Type>& d2){
	std::swap(d1.x_,d2.x_);
	std::swap(d1.dx_,d2.dx_);
	std::swap(d1.N_,d2.N_);
	std::swap(d1.conv_,d2.conv_);
	std::swap(d1.binning_,d2.binning_);
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
	d.set(x,dx,N,conv);
	return flux;
}

template<typename Type>
void Data<Type>::header_rst(std::string const& s, RST& rst) const {
	if(conv_){ rst.def(s,my::tostring(get_x())+" ("+my::tostring(get_dx())+")"); }
	else{ rst.def(s,"nc:"+my::tostring(get_x())+" ("+my::tostring(get_dx())+")"); }
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, Data<Type> const& d){
	if(w.is_binary()){ d.write(w); }
	else { w.stream()<<d; }
	return w;
}

template<typename Type>
IOFiles& operator>>(IOFiles& r, Data<Type>& d){
	if(r.is_binary()){ d = std::move(Data<Type>(r)); }
	else { r.stream()>>d; }
	return r;
}

template<typename Type>
void Data<Type>::write(IOFiles& w, std::string const& name) const {
	if(name != ""){
		w.add_header()->def(name,(conv_?"":"nc:")+my::tostring(get_x())+" ("+my::tostring(get_dx())+")"); 
	}
	w<<get_x()<<dx_<<N_<<conv_; 
	if(binning_){ 
		w<<true; 
		binning_->write(w);
	} else { w<<false; }
}
/*}*/

/*operator*/
/*{*/
template<typename Type>
Data<Type>& Data<Type>::operator=(Data<Type> d){
	swap_to_assign(*this,d);
	return (*this);
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void Data<Type>::add_sample(Data<Type> const& d){
	if(d.conv_){
		x_ += d.x_;
		dx_+= d.dx_;
		N_ += d.N_;
		conv_ = true;
	}
}

template<typename Type>
void Data<Type>::add_sample(){
	if(binning_){ binning_->add_sample(x_); }
	else { std::cerr<<"void Data<Type>::add_sample() : no binning"<<std::endl; }
}

template<typename Type>
void Data<Type>::merge(Data const& d){
	if(binning_ && d.binning_){
		binning_->merge(*d.binning_);
		N_++;
	} else { std::cerr<<"void Data<Type>::merge(Data const& d) : no binning"<<std::endl; }
}

template<typename Type>
void Data<Type>::compute_convergence(double const& tol) {
	if(binning_){ binning_->compute_convergence(tol,dx_,conv_);}
	else { std::cerr<<"void Data<Type>::compute_convergence(double const& tol) : no binning"<<std::endl; }
}

template<typename Type>
void Data<Type>::complete_analysis(double const& tol){
	if(binning_){ binning_->complete_analysis(tol,x_,dx_,conv_); }
	else { std::cerr<<"void Data<Type>::complete_analysis(double const& tol) : no Binning"<<std::endl; }
}

template<typename Type>
void Data<Type>::complete_analysis(){
	x_ = x_/N_;
	dx_ = dx_/(N_*sqrt(N_));
}

template<typename Type>
void Data<Type>::delete_binning(){
	if(binning_){
		delete binning_;
		binning_ = NULL;
	}
}
/*}*/
/*}*/

/*DataSet*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
DataSet<Type>::DataSet():
	size_(0),
	ds_(NULL)
{}

template<typename Type>
DataSet<Type>::DataSet(DataSet<Type> const& ds):
	size_(ds.size_),
	ds_(size_?new Data<Type>[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){
		ds_[i] = ds.ds_[i];
	}
	std::cout<<"dataset copy"<<std::endl;
}

template<typename Type>
DataSet<Type>::DataSet(DataSet<Type>&& ds):
	size_(ds.size_),
	ds_(ds.ds_)
{
	ds.ds_ = NULL;
}

template<typename Type>
DataSet<Type>::DataSet(IOFiles& r):
	size_(r.read<unsigned int>()),
	ds_(size_?new Data<Type>[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){
		ds_[i] = std::move(Data<Type>(r));
	}
}

template<typename Type>
void DataSet<Type>::set(unsigned int const& N){
	if(ds_){ delete[] ds_; }
	ds_ = new Data<Type>[N];
	size_ = N;
}

template<typename Type>
void DataSet<Type>::set(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].set(); }
}

template<typename Type>
void DataSet<Type>::set(unsigned int const& N, unsigned int const& B, unsigned int const& b, bool const& conv){
	set(N);
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
	rst.def(s,my::tostring(size_));
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
void DataSet<Type>::merge(DataSet<Type> const& ds){
	for(unsigned int i(0);i<size_;i++){ ds_[i].merge(ds[i]); }
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

template<typename Type>
void DataSet<Type>::delete_binning(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].delete_binning(); }
}
/*}*/
/*}*/
#endif
