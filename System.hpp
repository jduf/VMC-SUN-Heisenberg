#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "CreateSystem.hpp"
#include "Binning.hpp"
#include "Rand.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class System{
	public:
		/*!Create a System and extract the parameters from CreateSystem*/
		//{Description
		/*! Creates the system in function of the input parameters.
		 *
		 * - for each thread the system is independantly initialized
		 * - sets N, m, n, links_ and set tmp to the correct size 
		 * - allocates memory Ainv_
		 * - initialize the random number generator
		 */ //}
		System(CreateSystem* CS, unsigned int const& thread, unsigned int const& type);
		/*!Delete all the variables dynamically allocated*/
		virtual ~System();

		/*!Exchanges two particles of different color */
		virtual void swap();
		/*!Exchanges particle on site s1 with the one on site s2*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);
		/*!Virtual method that is called by MonteCarlo */
		virtual Type ratio()=0;
		/*!Updates only s_*/
		virtual void update();

		/*!Returns the status*/
		bool ready() const {return ready_;}
		/*!Returns the energy*/
		Data<double> get_energy() const {return E_.get_data();}
		/*!Returns the correlation*/
		Matrix<double> get_corr() const;
		/*!Returns the correlation*/
		Matrix<double> get_long_range_corr() const;

		//{Description
		/*!Computes the matrix element <a|H|b> where |a> and |b> differs by one
		 * permutation */
		//}
		void measure_new_step();	
		void add_sample();
		bool is_converged();
		
		void save(Write& w) const;

		void set();

		void operator>>(Write &w) const;

		/*!Pure virtual function that provides a way to check the System */
		virtual void print()=0;

	protected:
		unsigned int new_c[2];	//!< colors of the exchanged sites
		unsigned int new_s[2];	//!< sites that are exchanged
		unsigned int new_p[2];	//!< sites that are exchanged

		bool ready_;	

		unsigned int const N_;//!< colors' number
		unsigned int const n_;//!< sites' number
		unsigned int const m_;//!< particles per site' number
		unsigned int const M_;//!< particles' number of each color

		Matrix<unsigned int> s_;			//!< on the site i : s(i,0)=color, s(i,1)=row
		Matrix<unsigned int> const links_;	//!< list of links

		Rand* rnd_;	//!< generator of random numbers 

	private:
		Binning E_;
		Binning* corr_;				//!< correlation for each link 
		Binning* long_range_corr_;	//!< correlation for each link 

		/*!Check only if the new state has not the same color on one site*/
		bool is_new_state_forbidden();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
System<Type>::System(CreateSystem* CS, unsigned int const& thread, unsigned int const& type):
	ready_(false),
	N_(CS->get_N()),
	n_(CS->get_n()),
	m_(CS->get_m()),
	M_((m_*n_)/N_),
	s_(n_,m_),
	links_(CS->get_links()),
	rnd_(new Rand(100,thread))
{
	corr_ = new Binning[n_];
	if(type == 2){long_range_corr_ = new Binning[n_/3];}
	set();
	//E_.log("E");
	//long_range_corr_[n_/3-1].log("long_range_corr");
}

template<typename Type>
void System<Type>::set(){
	E_.set();
	for(unsigned int i(0);i<n_;i++){ corr_[i].set(); }
	if(long_range_corr_){
		for(unsigned int i(0);i<n_/3;i++){ long_range_corr_[i].set(); }
	}
}

template<typename Type>
System<Type>::~System(){
	//E_.plot();
	//long_range_corr_[n_/3-1].plot();

	if(rnd_){delete rnd_;}
	if(corr_){delete[] corr_;}
	if(long_range_corr_){delete[] long_range_corr_;}
}
/*}*/

/*public method*/
/*{*/
template<typename Type>
void System<Type>::swap(){
	new_s[0] = rnd_->get(n_);
	new_p[0] = rnd_->get(m_);
	new_c[0] = s_(new_s[0],new_p[0]);
	do {
		new_s[1] = rnd_->get(n_);
		new_p[1] = rnd_->get(m_);
		new_c[1] = s_(new_s[1],new_p[1]);
	} while(is_new_state_forbidden() || new_c[0] == new_c[1]);
}

template<typename Type>
void System<Type>::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	new_s[0] = s0;
	new_p[0] = p0;
	new_c[0] = s_(s0,p0);
	new_s[1] = s1;
	new_p[1] = p1;
	new_c[1] = s_(s1,p1);
}

template<typename Type>
void System<Type>::update(){
	/*update the sites*/
	s_(new_s[0],new_p[0]) = new_c[1];
	s_(new_s[1],new_p[1]) = new_c[0];
}

template<typename Type>
void System<Type>::measure_new_step(){
	E_ = 0.0;
	double r;
	for(unsigned int i(0);i<links_.row();i++){
		corr_[i] = 0.0;
		for(unsigned int p0(0); p0<m_; p0++){
			for(unsigned int p1(0); p1<m_; p1++){
				swap(links_(i,0),links_(i,1),p0,p1);
				if(!is_new_state_forbidden()){ 
					r = real(ratio());
					E_ += r; 
					corr_[i] += r;
				}
			}
		}
	}
	E_ /= n_;
	if(long_range_corr_){
		unsigned int x0(n_/3);
		for(unsigned int i(0);i<n_/3;i++){
			long_range_corr_[i] = 0.0;
			for(unsigned int p0(0); p0<m_; p0++){
				for(unsigned int p1(0); p1<m_; p1++){
					swap(x0,x0+i+1,p0,p1);
					if(!is_new_state_forbidden() && new_c[0]!=new_c[1]){ 
						long_range_corr_[i] += real(ratio());;
					}
				}
			}
		}
	}
}

template<typename Type>
void System<Type>::add_sample(){
	E_.add_sample();
	for(unsigned int i(0);i<n_;i++){
		corr_[i].add_sample();
	}
	if(long_range_corr_){
		for(unsigned int i(0);i<n_/3;i++){
			long_range_corr_[i].add_sample();
		}
	}
}

template<typename Type>
bool System<Type>::is_converged(){ 
	for(unsigned int i(0);i<n_;i++){
		corr_[i].is_converged(); 
	}
	if(long_range_corr_){
		for(unsigned int i(0);i<n_/3;i++){
			long_range_corr_[i].is_converged(); 
		}
	}
	return E_.is_converged(); 
}

template<typename Type>
void System<Type>::save(Write& w) const{
	//w("E (energy per site)",E_.get_mean());
	//w("DeltaE (absolute error)",E_.get_variance());
	w("corr (correlation on links)",get_corr());
	if(long_range_corr_){
		w("long_range_corr (long range correlations)",get_long_range_corr());
	}
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
bool System<Type>::is_new_state_forbidden(){
	for(unsigned int i(0); i<m_; i++){
		if(s_(new_s[0],i) == new_c[1] && i != new_p[0]){ return true; }
		if(s_(new_s[1],i) == new_c[0] && i != new_p[1]){ return true; }
	}
	return false;
}
/*}*/

template<typename Type>
Matrix<double> System<Type>::get_corr() const {
	Matrix<double> corr(n_,2);
	for(unsigned int i(0);i<n_;i++){
		//corr(i,0) = corr_[i].get_mean();
		//corr(i,1) = corr_[i].get_variance();
	}
	return corr;
}

template<typename Type>
Matrix<double> System<Type>::get_long_range_corr() const{
	Matrix<double> long_range_corr(n_/3,2);
	for(unsigned int i(0);i<n_/3;i++){
		//long_range_corr(i,0) = long_range_corr_[i].get_mean();
		//long_range_corr(i,1) = long_range_corr_[i].get_variance();
	}
	return long_range_corr;
}

template<typename Type>
void System<Type>::operator>>(Write &w) const{
	E_>>w;
	for(unsigned int i(0);i<n_;i++){
		corr_[i]>>w;
	}
	if(long_range_corr_){
		for(unsigned int i(0);i<n_/3;i++){
			long_range_corr_[i]>>w;
		}
	}
}
#endif
