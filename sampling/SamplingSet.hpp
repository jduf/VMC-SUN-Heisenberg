#ifndef DEF_SAMPLINGSET
#define DEF_SAMPLINGSET

#include "Sampling.hpp"

template<typename Type>
class CorrelatedSamplesSet;

template<typename Type>
class DataSet{
	public:
		/*!Default constructor*/
		DataSet();
		//void set(unsigned int const& N);
		~DataSet();

		void set(unsigned int const& N, bool const& conv);

		unsigned int size() const { return N_;}
		
		void add_sample(CorrelatedSamplesSet<Type> const& ss);
		void complete_analysis();

		Data<Type> const& operator[](unsigned int const& i) const 
		{assert(i<N_); return ds_[i];} 
		Data<Type>& operator[](unsigned int const& i) 
		{assert(i<N_); return ds_[i];} 

		void header_rst(std::string const& s, RST& rst) const;

	private:
		Data<Type>* ds_;
		unsigned int N_;
};

template<typename Type>
class CorrelatedSamplesSet{
	public:
		/*!Default constructor (needed if unused variable)*/
		CorrelatedSamplesSet(); 
		/*!Destructor*/
		~CorrelatedSamplesSet();
		/*!Set the class*/
		void set(unsigned int const& N, unsigned int const& B, unsigned int const& b);

		/*!Add sample to the bins*/
		void add_sample();
		/*!Run CorrelatedSamples::compute_convergence(double const& tol) on all set*/
		void compute_convergence(double const& tol);
		/*!Run CorrelatedSamples::complete_analysis(double const& tol) on all set*/
		void complete_analysis(double const& tol);

		CorrelatedSamples<Type> const& operator[](unsigned int const& i) const 
		{assert(i<N_); return cs_[i];}
		CorrelatedSamples<Type>& operator[](unsigned int const& i)
		{assert(i<N_); return cs_[i];}

		unsigned int size() const {return N_;}
		CorrelatedSamples<Type>* ptr() const {return cs_;}

		void header_rst(std::string const& s, RST& rst) const;

	private:
		CorrelatedSamples<Type>* cs_;
		unsigned int N_;
};

/*DataSet*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
DataSet<Type>::DataSet():
	ds_(NULL),
	N_(0)
{}

template<typename Type>
void DataSet<Type>::set(unsigned int const& N, bool const& conv){
	if(ds_){ delete[] ds_; }
	ds_ = new Data<Type>[N];
	N_ = N;
	for(unsigned int i(0);i<N_;i++){ ds_[i].set_conv(conv); }
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
	rst.def(s,tostring(N_));
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
		if(size != ds.size()){ ds.set(size,true); } 
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
void DataSet<Type>::add_sample(CorrelatedSamplesSet<Type> const& css){
	assert(N_ == css.size());
	for(unsigned int i(0);i<N_;i++){ ds_[i].add_sample(css[i]); }
}

template<typename Type>
void DataSet<Type>::complete_analysis(){
	for(unsigned int i(0);i<N_;i++){ ds_[i].complete_analysis(); }
}
/*}*/
/*}*/

/*CorrelationSamplesSet*/
/*{*/
/*constructors and destructor*/
/*{*/
template<typename Type>
CorrelatedSamplesSet<Type>::CorrelatedSamplesSet():
	cs_(NULL),
	N_(0)
{}

template<typename Type>
void CorrelatedSamplesSet<Type>::set(unsigned int const& N, unsigned int const& B, unsigned int const& b){
	if(cs_){ delete[] cs_; }
	cs_ = new CorrelatedSamples<Type>[N];
	N_ = N;
	for(unsigned int i(0);i<N_;i++){ cs_[i].set(B,b,false); }
}

template<typename Type>
CorrelatedSamplesSet<Type>::~CorrelatedSamplesSet(){
	if(cs_){delete[] cs_;}
}
/*}*/

/*i/o methods*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, CorrelatedSamplesSet<Type> const& css){
	for(unsigned int i(0);i<css.size();i++){ flux<<css[i]<<std::endl; }
	return flux;
}

template<typename Type>
void CorrelatedSamplesSet<Type>::header_rst(std::string const& s, RST& rst) const {
	rst.def(s,tostring(N_));
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, CorrelatedSamplesSet<Type> const& css){
	if(w.is_binary()){
		w<<css.size();
		for(unsigned int i(0);i<css.size();i++){ w<<css[i]; }
	} else {
		for(unsigned int i(0);i<css.size();i++){ w.stream()<<css[i]<<std::endl; }
	}
	return w;
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void CorrelatedSamplesSet<Type>::add_sample(){
	for(unsigned int i(0);i<N_;i++){ cs_[i].add_sample(); }
}

template<typename Type>
void CorrelatedSamplesSet<Type>::compute_convergence(double const& tol){
	for(unsigned int i(0);i<N_;i++){ cs_[i].compute_convergence(tol); }
}

template<typename Type>
void CorrelatedSamplesSet<Type>::complete_analysis(double const& tol){
	for(unsigned int i(0);i<N_;i++){ cs_[i].complete_analysis(tol); }
}
/*}*/
/*}*/
#endif
