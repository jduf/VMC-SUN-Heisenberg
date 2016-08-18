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
class Data{
	public:
		/*!Default constructor*/
		Data() = default;
		/*!Copy constructor*/
		Data<Type>(Data<Type> const& d);
		/*!Move constructor*/
		Data<Type>(Data<Type>&& d);
		/*!Constructor that reads from file*/
		Data<Type>(IOFiles& r);
		/*!Destructor*/
		~Data();
		void set(Type const& x, Type const& dx, double const& N, bool const& conv);
		void set(unsigned int const& B, unsigned int const& b, bool const& conv);
		void set();
		Data<Type>& operator=(Data<Type> d);

		void add_sample();
		void merge(Data const& d);

		void complete_analysis(double const& convergence_criterion);
		void delete_binning();

		double get_N() const { return binning_?binning_->get_N():N_; }
		Type get_x() const { return binning_?binning_->get_x():x_; }
		Type const& get_dx() const { return dx_; }
		bool const& get_conv() const { return conv_; }

		void set_x(Type const& x){ x_ = x; }
		void add(Type const& x){ x_ += x; }
		void divide(Type const& x){ x_ /= x; }

		std::string header_def() const { return (conv_?"":"nc:")+my::tostring(get_x())+" ("+my::tostring(get_dx())+")"; }
		void write(IOFiles& w) const;

	private:
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
				Binning() = delete;
				Binning& operator=(Binning) = delete;
				/*}*/

				/*!Set the class*/
				void set();

				/*!Add sample to the bins*/
				void add_sample(Type const& x);
				/*!Set x_ to the mean value, dx_ to the variance*/
				void complete_analysis(double const& convergence_criterion, Type& x, Type& dx, double& N, bool& conv);
				/*!Compute the mean value*/
				Type get_x() const { return ((Ml_(0)*m_bin_(0)*DPL_)+bin_[0](Ml_(0)))/get_N(); }
				/*!Compute the number of samples*/
				double get_N() const { return Ml_(0)*1.0*DPL_+dpl_; }
				/*!Merge this with b*/
				void merge(Binning const& b);

				void write(IOFiles& w) const;

			private:
				unsigned int const B_;//!< minimum number of biggest bins needed to compute variance
				unsigned int const b_;//!< l_+b_ rank of the biggest bin (b !> 30)
				unsigned int l_;	  //!< rank of the "smallest" bin
				unsigned int DPL_;	  //!< 2^l_ maximum number of element in the smallest bin l_
				unsigned int dpl_;	  //!< current number of element in each bin of rank l_

				Vector<unsigned int> Ml_;//!< number filled (or partially filled) bins of rank l : Ml_(l) = Ml_(0)/2^l
				Vector<Type> m_bin_;	 //!< mean of the Binnings
				Vector<Type>*  bin_;	 //!< binnings

				bool recompute_dx_usefull_;	//!< true if dx should be recomputed

				/*!Recursive method that add samples in the different bins*/
				void add_bin(unsigned int l, Type const& a, Type const& b);
				void do_merge(Vector<Type> const& bin, unsigned int const& dpl, unsigned int const& DPL, unsigned int const& Ml);
		};

		Type x_			  = 0.0;
		Type dx_		  = 0.0;
		double N_		  = 0.0;
		bool conv_		  = false;
		Binning* binning_ = NULL;

		void swap_to_assign(Data<Type>& d1, Data<Type>& d2);
};

template<typename Type>
class DataSet{
	public:
		/*!Default constructor*/
		DataSet() = default;
		/*!Copy constructor*/
		DataSet<Type>(DataSet<Type> const& d);
		/*!Copy constructor*/
		DataSet<Type>(DataSet<Type>&& d);
		/*!Constructor that reads from file*/
		DataSet<Type>(IOFiles& r);
		/*!Destructor*/
		~DataSet();
		void set(unsigned int const& size, unsigned int const& B, unsigned int const& b, bool const& conv);
		void set(unsigned int const& size);
		void set();
		DataSet<Type>& operator=(DataSet<Type> d);

		void add_sample();
		void merge(DataSet<Type> const& ds);

		void complete_analysis(double const& convergence_criterion);
		void delete_binning();

		unsigned int size() const { return size_; }

		Data<Type> const& operator[](unsigned int const& i) const
		{assert(i<size_); return ds_[i]; }
		Data<Type>& operator[](unsigned int const& i)
		{assert(i<size_); return ds_[i]; }

		std::string header_def() const { return "DataSet("+my::tostring(size_)+")"; }

	private:
		unsigned int size_ = 0;
		Data<Type>* ds_    = NULL;

		void swap_to_assign(DataSet<Type>& ds1, DataSet<Type>& ds2);
};

/*Binning*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Data<Type>::Binning::Binning(unsigned int const& B, unsigned int const& b):
	B_(B),
	b_(b),
	bin_(new Vector<Type>[b])
{ set(); }

template<typename Type>
Data<Type>::Binning::Binning(Binning const& b):
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
}

template<typename Type>
Data<Type>::Binning::Binning(Binning&& b):
	B_(b.B_),
	b_(b.b_),
	l_(b.l_),
	DPL_(b.DPL_),
	dpl_(b.dpl_),
	Ml_(b.Ml_),
	m_bin_(std::move(b.m_bin_)),
	bin_(b.bin_),
	recompute_dx_usefull_(b.recompute_dx_usefull_)
{ b.bin_ = NULL; }

template<typename Type>
Data<Type>::Binning::Binning(IOFiles& r):
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
void Data<Type>::Binning::set(){
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
Data<Type>::Binning::~Binning(){
	if(bin_){ delete[] bin_; }
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void Data<Type>::Binning::add_sample(Type const& x){
	/*!add new entry to the bins*/
	if(DPL_ == ++dpl_){
		add_bin(0,2.0/DPL_*x,2.0/DPL_*bin_[0](Ml_(0)));
		recompute_dx_usefull_ = true;
		dpl_ = 0;
		/*!update the bins if the bigger bin is big enough*/
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
void Data<Type>::Binning::merge(Binning const& other){
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
			dpl_= other.dpl_;
			DPL_= other.DPL_;
			l_  = other.l_;
			m_bin_ = other.m_bin_;
			recompute_dx_usefull_ = other.recompute_dx_usefull_;
			do_merge(tmp_bin,tmp_dpl,tmp_DPL,tmp_Ml);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : B_ != b.B_ || b_ != b.b_ "<<std::endl; }
}

template<typename Type>
void Data<Type>::Binning::complete_analysis(double const& convergence_criterion, Type& x, Type& dx, double& N, bool& conv){
	x = get_x();
	N = get_N();
	if(l_>30){ std::cerr<<__PRETTY_FUNCTION__<<" : possibility of 'unsigned int' overflow : l_="<<l_<<std::endl; }
	if(recompute_dx_usefull_ && Ml_(b_-1)>=B_){
		recompute_dx_usefull_ = false;
		/*!compute the variance for each bin*/
		Vector<Type> var_bin(b_,0.0);
		for(unsigned int l(0);l<b_;l++){
			for(unsigned int j(0);j<Ml_(l);j++){
				var_bin(l) += my::norm_squared(bin_[l](j)-m_bin_(l));
			}
			var_bin(l) = sqrt(var_bin(l) / (Ml_(l)*(Ml_(l)-1)));
		}

		/*!do a linear regression and if the slope is almost flat, the system
		 * is believed to be converged*/
		Vector<Type> vx(b_);
		Type xb(0); //mean of x
		Type yb(0); //mean of y
		Type num(0);
		Type den(0);
		for(unsigned int i(0);i<b_;i++){
			vx(i) = l_+i;
			xb += vx(i);
			yb += var_bin(i);
		}
		xb /= b_;
		yb /= b_;
		for(unsigned int i(0);i<b_;i++){
			num += (vx(i)-xb)*(var_bin(i)-yb);
			den += (vx(i)-xb)*(vx(i)-xb);
		}

		/*!set the incertitude over x as the mean over all bins of the variance
		 * if the slope is positive increase the incertitude */
		dx = yb;
		if(num/den>0.0){ dx += num/den * xb; }

		/*!if the slope is flat enough, the convergence is reached*/
		if(std::abs(num/(den*x))<convergence_criterion){ conv = true; }
		else{ conv = false; }
	}
}

template<typename Type>
void Data<Type>::Binning::write(IOFiles& w) const {
	w<<B_<<b_<<l_<<DPL_<<dpl_<<Ml_<<m_bin_<<bin_[0];
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void Data<Type>::Binning::add_bin(unsigned int l, Type const& a, Type const& b){
	bin_[l](Ml_(l)) = (a+b)/2.0;
	m_bin_(l) = (m_bin_(l)*double(Ml_(l))+bin_[l](Ml_(l)))/(Ml_(l)+1.0);
	Ml_(l)++;
	/*!Create the next (bigger) bin*/
	if(Ml_(l)%2==0 && l<b_-1){
		add_bin(l+1,bin_[l](Ml_(l)-1),bin_[l](Ml_(l)-2));
	}
}

template<typename Type>
void Data<Type>::Binning::do_merge(Vector<Type> const& bin, unsigned int const& dpl, unsigned int const& DPL, unsigned int const& Ml){
	Type old_last_bin(bin_[0](Ml_(0)));
	unsigned int old_dpl(dpl_);
	bin_[0](Ml_(0)) = 0;
	dpl_ = 0;
	for(unsigned int i(0);i<Ml;i++){
		dpl_ += DPL-1;
		add_sample(DPL*bin(i));
	}
	Type tmp((bin(Ml)+old_last_bin)/(dpl+old_dpl));
	for(unsigned int i(0);i<dpl+old_dpl;i++){ add_sample(tmp); }
}
/*}*/
/*}*/

/*Data*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Data<Type>::Data(Data<Type> const& d):
	x_(d.x_),
	dx_(d.dx_),
	N_(d.N_),
	conv_(d.conv_),
	binning_(d.binning_?new Data<Type>::Binning(*d.binning_):NULL)
{}

template<typename Type>
Data<Type>::Data(Data<Type>&& d):
	x_(d.x_),
	dx_(d.dx_),
	N_(d.N_),
	conv_(d.conv_),
	binning_(d.binning_)
{ d.binning_ = NULL; }

template<typename Type>
Data<Type>::Data(IOFiles& r):
	x_(r.read<Type>()),
	dx_(r.read<Type>()),
	N_(r.read<double>()),
	conv_(r.read<bool>()),
	binning_(r.read<bool>()?new Data<Type>::Binning(r):NULL)
{}

template<typename Type>
void Data<Type>::set(Type const& x, Type const& dx, double const& N, bool const& conv){
	x_ = x;
	dx_ = dx;
	N_ = N;
	conv_ = conv;
}

template<typename Type>
void Data<Type>::set(unsigned int const& B, unsigned int const& b, bool const& conv){
	if(binning_){ delete binning_; }
	binning_ = new Data<Type>::Binning(B,b);
	conv_ = conv;
	x_ = 0.0;
	dx_ = 0.0;
	N_ = 0.0;
}

template<typename Type>
void Data<Type>::set(){
	binning_->set();
	x_ = 0.0;
	dx_ = 0.0;
	N_ = 0.0;
	conv_ = false;
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
	double N(0.0);
	bool conv(false);
	flux>>x>>dx>>N>>conv;
	d.set(x,dx,N,conv);
	return flux;
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
void Data<Type>::write(IOFiles& w) const {
	/*!the call of get_x() and get_N() is useful because x_ and N_ are not
	 * necessarily set to their correct value*/
	w<<get_x()<<dx_<<get_N()<<conv_;
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
void Data<Type>::add_sample(){
	if(binning_){ binning_->add_sample(x_); }
	else { std::cerr<<__PRETTY_FUNCTION__<<" : no binning"<<std::endl; }
}

template<typename Type>
void Data<Type>::merge(Data const& d){
	if(d.binning_){
		if(!binning_){ binning_ = new Data<Type>::Binning(*d.binning_); }
		else { binning_->merge(*d.binning_); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no binning"<<std::endl; }
}

template<typename Type>
void Data<Type>::complete_analysis(double const& convergence_criterion){
	if(binning_){ binning_->complete_analysis(convergence_criterion,x_,dx_,N_,conv_); }
	else { std::cerr<<__PRETTY_FUNCTION__<<" : no binning"<<std::endl; }
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
DataSet<Type>::DataSet(DataSet<Type> const& ds):
	size_(ds.size_),
	ds_(size_?new Data<Type>[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){
		ds_[i] = ds.ds_[i];
	}
}

template<typename Type>
DataSet<Type>::DataSet(DataSet<Type>&& ds):
	size_(ds.size_),
	ds_(ds.ds_)
{ ds.ds_ = NULL;
	//std::cerr<<"move DataSet"<<std::endl;
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
void DataSet<Type>::set(unsigned int const& size){
	if(ds_){ delete[] ds_; }
	ds_ = new Data<Type>[size];
	size_ = size;
}

template<typename Type>
void DataSet<Type>::set(unsigned int const& size, unsigned int const& B, unsigned int const& b, bool const& conv){
	set(size);
	for(unsigned int i(0);i<size_;i++){ ds_[i].set(B,b,conv); }
}

template<typename Type>
void DataSet<Type>::set(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].set(); }
}

template<typename Type>
DataSet<Type>::~DataSet(){
	if(ds_){ delete[] ds_; }
}

template<typename Type>
void DataSet<Type>::swap_to_assign(DataSet<Type>& ds1, DataSet<Type>& ds2){
	std::swap(ds1.size_,ds2.size_);
	std::swap(ds1.ds_,ds2.ds_);
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
	std::cerr<<__PRETTY_FUNCTION__<<" :  not defined"<<std::endl;
	return flux;
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

/*operator*/
/*{*/
template<typename Type>
DataSet<Type>& DataSet<Type>::operator=(DataSet<Type> ds){
	swap_to_assign(*this,ds);
	return (*this);
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void DataSet<Type>::add_sample(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].add_sample(); }
}

template<typename Type>
void DataSet<Type>::merge(DataSet<Type> const& ds){
	if(size_ != ds.size_){ set(ds.size_); }
	for(unsigned int i(0);i<size_;i++){ ds_[i].merge(ds[i]); }
}

template<typename Type>
void DataSet<Type>::complete_analysis(double const& convergence_criterion){
	for(unsigned int i(0);i<size_;i++){ ds_[i].complete_analysis(convergence_criterion); }
}

template<typename Type>
void DataSet<Type>::delete_binning(){
	for(unsigned int i(0);i<size_;i++){ ds_[i].delete_binning(); }
}
/*}*/
/*}*/
#endif
