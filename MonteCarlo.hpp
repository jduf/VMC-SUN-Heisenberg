#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"
#include "Write.hpp"

//{Description
/*! Class MonteCarlo
 *
 * Implement the Monte-Carlo algorithm. This alorithm let a system (System*)
 * evolves according to a Markov process. To work properly, System* has to be
 * written with :
 * 
 * + System->ratio() : compute the probability to accept the next configuration
 * + System->update() : update the old cufiguration to the new one
 * + System->measure() : measure an observable according to the current configuration
 *
 * Each time that a class is instanciated, a random number generator is created
 * according to the thread on which the code is running. The same thread number
 * is then transmitted to create the (System*). Once the (*System) is created,
 * it is thermilized.
 *
 * The run(what,N_MC) method lunches the Monte-Carlo simulation.
 */
//}
template <typename Type>
class MonteCarlo{
	public:
		/*!The container must contain a filename to store the executing data,
		 * the parameter required for (*System) and an optional t_max option*/
		MonteCarlo(Container const& input); 
		~MonteCarlo();

		/*!Run the Monte-Carlo algorithme */
		void run(unsigned int what, unsigned int N_MC=0);
		/*!Saves the essential data in the "result" file*/
		void save(Container& result);

		/*!Get the energy per site*/
		double get_energy() const { return this->E_; };
		unsigned int get_status() const { return status; };

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc); 
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		/*!Find the next configuration and measure it*/
		void next_step(double& E_step, Vector<double>& corr_step);
		//{Private methods that give a shutoff condition
		/*!Stops the simulation when
		 *
		 * - convergence is reached
		 * - time limit is up
		 * - if the energy diverges
		 *
		 * ands save each step in the "output" file
		 */
		//}
		bool keepon(double E_step, unsigned int N_MC);
		void add_bin(unsigned int l, double a, double b);
		void compute_error();

		System<Type>* S; 		//!< Pointer to a Fermionic or Bosonic System 
		Rand* rnd; 				//!< Pointer to a random number generator
		Time timer; 			//!< To stop the simulation after time_limit seconds
		unsigned int t_max;		//!< Time limit in second, by default 5min
		unsigned int status; 	//!< Not Lunched:0 Lunched:1 Successful:2 Time elapsed:3

		unsigned int N_step_;	//!< Number of steps already done
		unsigned int B_;		//!< Minimum number of biggest bins needed to compute variance
		unsigned int b_;		//!< Number of different binning
		unsigned int dpl_;		//!< 2^l, l:number of element in the smallest bin
		unsigned int dpb_;

		Vector<unsigned int> Ml_;//!< Number of element in binning of rank l
		Vector<double>*bin_;	//!< Binnings
		Vector<double> m_bin_;	//!< Mean of the Binnings
		Vector<double> usl_;	//!< 

		Vector<double> corr_;
		double E_; 				//!< Value that the MC algorithm tries to compute
		double DeltaE_;			//!< Error on E_

		std::string filename_;	//!< Filename for the output
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(Container const& input):
	S(NULL),
	rnd(NULL),
	t_max(5*60),
	status(0),
	N_step_(0),
	B_(50),
	b_(5),
	dpl_(2),
	dpb_(1),
	Ml_(b_,0),
	bin_(new Vector<double>[b_]),
	m_bin_(b_,0.0),
	usl_(b_,0.0),
	E_(0),
	DeltaE_(0),
	filename_("")
{
	for(unsigned int i(b_);i>0;i--){
		dpb_ *= 2;
		bin_[i-1].set(B_*dpb_,0);
	}
	for(unsigned int i(0);i<b_;i++){
		usl_(i) = 1.0/(i+1);
	}
	if(input.check("filename")){
		input.get("filename",filename_);
	}
	unsigned int thread(omp_get_thread_num());
	rnd=new Rand(1e4,thread);
	if(input.get<Vector<unsigned int> >("ref")(0)){ S=new SystemFermionic<Type>; }
	else { S=new SystemBosonic<Type>; }
	input.get("t_max",t_max);
	status = S->init(input,thread);
	if(status){
		double ratio(0.0);
		for(unsigned int i(0);i<1e5;i++){
			S->swap();
			ratio = norm_squared(S->ratio());
			if( ratio > rnd->get() ){
				S->update();
			}
		}
		corr_.set(S->sts_.row());
	}
}

template<typename Type>
MonteCarlo<Type>::~MonteCarlo(){
	delete S;
	delete rnd;
	delete bin_;
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::run(unsigned int what, unsigned int N_MC){
	if(status){
		double E_step(0.0);
		Vector<double> corr_step(S->sts_.row());
		S->measure(E_step,corr_step);
		switch(what){
			case 1: 
				{
					//Write output(filename_+".out");
					do{
						next_step(E_step,corr_step); 
						//output<<N_step
							//<<" "<<E_ 
							//<<" "<<err
							//<<Write::endl;
					} while(keepon(E_step,N_MC));
					//Write correlation_file(filename_+".corr");
					//correlation_file<<corr_<<Write::endl;
				} break;
			case 2:
				{
					do{ next_step(E_step,corr_step); }
					while(keepon(E_step,N_MC));
				} break;
		}
		E_ /= N_step_;
	}
}

template<typename Type>
void MonteCarlo<Type>::save(Container& result){
	result.set("N_step",N_step_);
	result.set("E",E_);
	result.set("dE",DeltaE_);
	result.set("status",status);
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::next_step(double& E_step, Vector<double>& corr_step){
	S->swap();
	if( norm_squared(S->ratio()) > rnd->get() ){
		S->update();
		S->measure(E_step,corr_step);
	}
	N_step_++;
	E_ += E_step;
	corr_ += corr_step;
}

template<typename Type>
bool MonteCarlo<Type>::keepon(double E_step, unsigned int N_MC){
	if(b_>30){std::cerr<<"will bug"<<std::endl;}
	/*!add new entry to the bins*/
	if(N_step_%dpl_==0){//l_=1 at the very beginning
		add_bin(0,E_step,bin_[0](Ml_(0))/(dpl_-1));
	} else {
		bin_[0](Ml_(0)) += E_step;
	}

	/*!update the bins if the bigger binning is big enough*/
	if(N_step_==B_*dpl_*dpb_){//B*2^(l+b)
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
	}

	/*!compute the variance and the running condition if one bin is added*/
	if(N_step_%dpl_==0 && Ml_(b_-1) >= B_ ){
		bool kpn(true);
		if( std::abs(E_)/N_step_ > 1e2 ){ 
			std::cerr<<"the simulation should be restarted after "<<N_step_<<" steps, E="<<E_<<std::endl;
			kpn=false;
		}
		if(timer.limit_reached(t_max)){
			std::cerr<<"time limit reached : stopping simulation"<<std::endl;
			kpn=false;
			status = 3;
		}
		if(N_step_ > N_MC && N_MC != 0){
			std::cerr<<"N_step>N_MC : stopping simulation"<<std::endl;
			kpn=false;
			status = 3;
		}
		compute_error();
		if(std::abs(DeltaE_/E_)*N_step_<1e-3){
			kpn=false;
			status = 2;
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
