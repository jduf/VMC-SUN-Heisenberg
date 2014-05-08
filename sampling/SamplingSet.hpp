#ifndef DEF_SAMPLINGSET
#define DEF_SAMPLINGSET

#include "Sampling.hpp"

template<typename Type>
class CorrelatedSamplesSet{
	public:
		/*!Default constructor (needed if unused variable)*/
		CorrelatedSamplesSet(); 
		/*!Destructor*/
		~CorrelatedSamplesSet();

		void set(unsigned int const& N, unsigned int const& B, unsigned int const& b);

		CorrelatedSamples<Type> const& operator[](unsigned int const& i) const 
		{assert(i<N_); return cs_[i];}
		CorrelatedSamples<Type>& operator[](unsigned int const& i)
		{assert(i<N_); return cs_[i];}

		void compute_convergence(double const& tol);
		void add_sample();

		unsigned int size() const {return N_;}

		CorrelatedSamples<Type>* ptr() const {return cs_;}

	private:
		CorrelatedSamples<Type>* cs_;
		unsigned int N_;
};

template<typename Type>
class DataSet{
	public:
		/*!Default constructor*/
		DataSet();
		/*!Constructor using set*/
		DataSet(unsigned int const& N);
		//void set(unsigned int const& N);
		~DataSet();

		void add_sample(CorrelatedSamplesSet<Type> const& ss);

		Data<Type> const& operator[](unsigned int const& i) const 
		{assert(i<N_); return ds_[i];} 
		Data<Type>& operator[](unsigned int const& i) 
		{assert(i<N_); return ds_[i];} 

		unsigned int size() const { return N_;}
		
	private:
		Data<Type>* ds_;
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
DataSet<Type>::DataSet(unsigned int const& N):
	ds_(new Data<Type>[N]),
	N_(N)
{}

//template<typename Type>
//void DataSet<Type>::set(unsigned int const& N){
	//ds_ = new Data<Type>[N];
	//N_ = N;
//}

template<typename Type>
DataSet<Type>::~DataSet(){
	if(ds_){ delete[] ds_;}
}
/*}*/

//template<typename Type>
//DataSet<Type>& DataSet<Type>::operator=(CorrelatedSamplesSet<Type> const& css){
	//if(ds_){ delete[] ds_;}
	//ds_ = new Data<Type>[css.size()];
	//N_ = css.size();
	//for(unsigned int i(0);i<N_;i++){ ds_[i] = css[i]; }
	//return (*this);
//}

//template<typename Type>
//DataSet<Type>& DataSet<Type>::operator+=(CorrelatedSamplesSet<Type> const& css){
	//assert(N_ == css.size());
	//for(unsigned int i(0);i<N_;i++){ ds_[i] += css[i]; }
	//return (*this);
//}

/*public methods that modify the class*/
/*{*/
template<typename Type>
void DataSet<Type>::add_sample(CorrelatedSamplesSet<Type> const& css){
	assert(N_ == css.size());
	for(unsigned int i(0);i<N_;i++){ ds_[i].add_sample(css[i]); }
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
	if(cs_){delete[] cs_;}
	cs_ = new CorrelatedSamples<Type>[N];
	N_ = N;
	for(unsigned int i(0);i<N_;i++){cs_[i].set(B,b);}
}

template<typename Type>
CorrelatedSamplesSet<Type>::~CorrelatedSamplesSet(){
	if(cs_){delete[] cs_;}
}
/*}*/

/*public methods that modify the class*/
/*{*/
template<typename Type>
void CorrelatedSamplesSet<Type>::add_sample(){
	for(unsigned int i(0);i<N_;i++){cs_[i].add_sample();}
}

template<typename Type>
void CorrelatedSamplesSet<Type>::compute_convergence(double const& tol){
	for(unsigned int i(0);i<N_;i++){cs_[i].compute_convergence(tol);}
}
/*}*/
/*}*/
#endif
