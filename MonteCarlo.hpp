#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"

/*{! Class MonteCarlo
 *
 * Implement the Monte-Carlo algorithm. This alorithm let a System evolves
 * according to a Markov process. To work properly, System has to contain at
 * least these methods : 
 *
 * - System::ratio() : compute the probability to accept the next configuration
 * - System::update() : update the old cufiguration to the new one
 * - System::measure() : measure an observable according to the current configuration 
 *
 * The MonteCarlo class contains a pointer to System. According to the system,
 * the constructor creates a SystemFermionic or SystemBosonic. As the
 * constructors of these classes need a pointer to a CreateSystem, the
 * constructor of MonteCarlo must provide a pointer on the same class.
 *
 * Each time that a class is instanciated, a random number generator is created
 * according to the thread on which the code is running. The same thread number
 * is then transmitted to create the System.
 *
 * Once the System is created, it is thermalized.
 *
 * The MonteCarlo::run() method lunches the Monte-Carlo simulation.
 }*/
template <typename Type>
class MonteCarlo{
	public:
		/*!Constructor, if Nmaxsteps_ isn't set, Nsteps_ may overflow*/
		MonteCarlo(CreateSystem* CS, unsigned int tmax, unsigned int Nmaxsteps=1e9); 
		/*!Simple destructor*/
		~MonteCarlo();

		/*!Run the Monte-Carlo algorithm*/
		void run();
		/*!Saves E_, DeltaE_, status_, Nsteps_ and corr_ in file w*/
		void save(Write& w) const;

		/*!Get the energy per site*/
		double get_energy() const { return E_; };
		/*!Get the abolute error per site*/
		double get_error() const { return DeltaE_; };
		/*!Get the correlation between neighbouring sites*/
		Vector<double> get_corr() const { return corr_; };
		/*!Get the long range correlation*/
		Vector<double> get_long_range_corr() const { return long_range_corr_; };
		/*!Get the status of the simulation*/
		unsigned int get_status() const { return status_; };

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc); 
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		/*!Initialize*/
		void init();
		/*!Find the next configuration and measure it*/
		void next_step(double& E_step, Vector<double>& corr_step, Vector<double>& long_range_corr_step);
		//{Private method that gives a shutoff condition
		/*!Stops the simulation when
		 *
		 * - convergence is reached
		 * - time limit is up
		 * - if the energy diverges
		 *
		 * If the E_ diverges, the simulation is restarted
		 */
		//}
		bool keepon(double E_step);
		/*!Recursive method that add datas in the different bins*/
		void add_bin(unsigned int l, double a, double b);
		/*!Compute the error according to Mathias Troyer*/
		void compute_error();

		System<Type>* S; 		//!< Pointer to a Fermionic or Bosonic System 
		Rand* rnd; 				//!< Pointer to a random number generator
		Time timer; 			//!< To stop the simulation after time_limit seconds
		unsigned int const tmax_;//!< Time limit in second, by default 5min
		unsigned int const Nmaxsteps_;//!< Maximum number of steps allowed
		unsigned int Nsteps_;	//!< Number of steps already done
		unsigned int status_; 	//!< Degenerate :0 No initial state:1 Launched:2 Time:3 Step:4 Successful:5 
		std::string status_message_;//!< Output message that explicits status_

		unsigned int const B_;	//!< Minimum number of biggest bins needed to compute variance
		unsigned int const b_;	//!< Number of different binning (b !> 30)
		unsigned int dpl_;		//!< 2^l, l:number of element in the smallest bin
		unsigned int B2pbl_;	//!< B2^(l+b)

		Vector<unsigned int> Ml_;//!< Number of element in binning of rank l
		Vector<double> usl_;	//!< 1/l, usefull for the linear regression
		Vector<double> m_bin_;	//!< Mean of the Binnings
		Vector<double>*bin_;	//!< Binnings

		double E_; 				//!< Value that the MC algorithm tries to compute
		double DeltaE_;			//!< Error on E_
		double const tol_;		//!< If DeltaE_/E_ < tol_, stops the simulation
		Vector<double> corr_;	//!< Correlation for each link 
		Vector<double> long_range_corr_;//!< Correlation for each link 
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(CreateSystem* CS, unsigned int tmax, unsigned int Nmaxsteps):
	S(NULL),
	rnd(NULL),
	tmax_(tmax),
	Nmaxsteps_(Nmaxsteps),
	status_(0),
	status_message_(""),
	B_(50),
	b_(5),
	bin_(new Vector<double>[b_]),
	tol_(1e-3),
	corr_(CS->get_num_links(),0),
	long_range_corr_(CS->get_n()/3-1,0)
{
	init();
	unsigned int thread(omp_get_thread_num());
	rnd=new Rand(1e4,thread);
	if(CS->is_bosonic()){ S=new SystemBosonic<Type>(CS,thread); }
	else { S=new SystemFermionic<Type>(CS,thread); }
	status_ = S->get_status();
	switch(status_){
		case 0:{status_message_ = "degenerate";}break;
		case 1:{status_message_ = "no initial state found";}break;
		case 2:{
				   double ratio(0.0);
				   for(unsigned int i(0);i<1e5;i++){
					   S->swap();
					   ratio = norm_squared(S->ratio());
					   if( ratio > rnd->get() ){
						   S->update();
					   }
				   }
				   status_message_ = "system thermalized";
			   }break;
		default:{std::cerr<<"MonteCarlo : System::status unknown"<<std::endl;}
	}
}

template<typename Type>
MonteCarlo<Type>::~MonteCarlo(){
	delete S;
	delete rnd;
	delete[] bin_;
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::run(){
	if(status_ == 2){/*passed the first two steps*/
		double E_step(0.0);
		Vector<double> corr_step(corr_);
		Vector<double> long_range_corr_step(long_range_corr_);
		S->measure(E_step,corr_step);
		S->long_range_corr(long_range_corr_step);
		do{ next_step(E_step,corr_step,long_range_corr_step); }
		while(keepon(E_step));
		corr_ /= Nsteps_;
		long_range_corr_ /= Nsteps_;
	}
}

template<typename Type>
void MonteCarlo<Type>::save(Write& w) const{
	w("E (energy per site)",E_);
	w("DeltaE (absolute error)",DeltaE_);
	w("Nsteps (number of steps)",Nsteps_);
	w("status ("+status_message_+")",status_);
	w("corr (correlation on links)",corr_);
	w("long_range_corr (correlation on links)",long_range_corr_);
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::init(){
	E_ = 0;
	DeltaE_ = 0;
	corr_.set(corr_.size(),0);
	Nsteps_ = 0;
	dpl_ = 2;
	B2pbl_ = B_;
	Ml_.set(b_,0);
	m_bin_.set(b_,0);
	usl_.set(b_,0.0);
	for(unsigned int i(b_);i>0;i--){
		B2pbl_ *= 2;
		bin_[i-1].set(B2pbl_,0);
	}
	for(unsigned int i(0);i<b_;i++){
		usl_(i) = 1.0/(i+1);
	}
}

template<typename Type>
void MonteCarlo<Type>::next_step(double& E_step, Vector<double>& corr_step, Vector<double>& long_range_corr_step){
	S->swap();
	if( norm_squared(S->ratio()) > rnd->get() ){
		S->update();
		S->measure(E_step,corr_step);
		S->long_range_corr(long_range_corr_step);
	}
	E_ *= Nsteps_;
	E_ += E_step;
	Nsteps_++;
	E_ /= Nsteps_;
	corr_ += corr_step;
	long_range_corr_ += long_range_corr_step;
}

template<typename Type>
bool MonteCarlo<Type>::keepon(double E_step){
	/*!add new entry to the bins*/
	if(Nsteps_%dpl_==0){//l_=1 at the very beginning
		add_bin(0,E_step,bin_[0](Ml_(0))/(dpl_-1));
	} else {
		bin_[0](Ml_(0)) += E_step;
	}

	/*!update the bins if the bigger binning is big enough*/
	if(Nsteps_==B2pbl_){//B*2^(l+b)
		for(unsigned int l(0);l<b_-1;l++){
			for(unsigned int j(0);j<bin_[l+1].size();j++){
				bin_[l](j) = bin_[l+1](j);
			}
			for(unsigned int j(bin_[l+1].size());j<bin_[l].size();j++){
				bin_[l](j) = 0.0;
			}
			Ml_(l) = Ml_(l+1);
			m_bin_(l) = m_bin_(l+1);
			usl_(l) = usl_(l+1);
		}
		Ml_(b_-1)=0;
		bin_[b_-1].set(2*B_,0);
		usl_(b_-1)=usl_(b_-1)/2.0;
		for(unsigned int j(0);j<B_;j++){
			add_bin(b_-1,bin_[b_-2](j+1),bin_[b_-2](j));
		}
		dpl_*=2;
		B2pbl_*=2;
	}

	/*!compute the variance and the running condition if one bin is added*/
	if(Nsteps_%dpl_==0 && Ml_(b_-1) >= B_ ){
		bool kpn(true);
		if( std::abs(E_) > 1e2 ){ 
			std::cerr<<"Simulation diverges => is restarted"<<std::endl;
			init();
		}
		if(timer.limit_reached(tmax_)){
			status_message_ = "time limit reached = " + tostring(tmax_);
			kpn=false;
			status_ = 3;
		}
		if(Nsteps_ > Nmaxsteps_){
			status_message_ = "Nsteps>Nmaxsteps = " + tostring(Nmaxsteps_);
			kpn=false;
			status_ = 4;
		}
		compute_error();
		if(std::abs(DeltaE_/E_)<tol_){
			status_message_ = "converged = " + tostring(DeltaE_/E_) + "<" + tostring(tol_);
			kpn=false;
			status_ = 5; /*3rd step successful*/
		}
		return kpn;
	} else {
		return true;
	}
}

template<typename Type>
void MonteCarlo<Type>::add_bin(unsigned int l, double a, double b){
	bin_[l](Ml_(l)) = (a+b)/2.0;
	m_bin_(l) = (m_bin_(l)*Ml_(l)+bin_[l](Ml_(l)))/(Ml_(l)+1);
	Ml_(l)++;
	/*!Create the next (bigger) bin*/
	if(Ml_(l)%2==0 && l<b_-1){
		add_bin(l+1,bin_[l](Ml_(l)-1),bin_[l](Ml_(l)-2));
	}
}

template<typename Type>
void MonteCarlo<Type>::compute_error(){
	/*!Compute the variance for each bin*/
	Vector<double> var_bin(b_,0.0);
	for(unsigned int l(0);l<b_;l++){
		for(unsigned int j(0);j<Ml_(l);j++){
			var_bin(l) += (bin_[l](j)-m_bin_(l))*(bin_[l](j)-m_bin_(l));
		}
		var_bin(l) = sqrt(var_bin(l) /(Ml_(l)*(Ml_(l)-1)));
	}

	/*!Do a linear regression to get the varience in the limit l->infty*/
	double xb(0);
	double yb(0);
	for(unsigned int i(0);i<b_;i++){
		xb += usl_(i);
		yb += var_bin(i);
	}
	xb /= b_;
	yb /= b_;
	double num(0);
	double den(0);
	for(unsigned int i(0);i<b_;i++){
		num += (usl_(i)-xb)*(var_bin(i)-yb);
		den += (usl_(i)-xb)*(usl_(i)-xb);
	}
	DeltaE_ = yb-num/den*xb;
}
/*}*/
#endif
