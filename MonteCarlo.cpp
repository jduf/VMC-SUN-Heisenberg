#include "MonteCarlo.hpp"

/*constructors and destructor*/
/*{*/
MonteCarlo::MonteCarlo(MCSystem* S, unsigned int const& tmax):
	tmax_(tmax),
	S_(S),
	rnd_(0.0,1.0)
{}
/*}*/

/*public methods*/
/*{*/
void MonteCarlo::thermalize(unsigned int const& thermalization_steps){
	if(!S_->get_status()){
		for(unsigned int i(0);i<thermalization_steps;i++){
			S_->swap();
			ratio_ = my::norm_squared(S_->ratio());
			if( ratio_ > 1.0 || ratio_ > rnd_.get() ){ S_->update(); }
		}
		S_->measure_new_step();
	}
}

void MonteCarlo::run(){
	time_.set();
	if(!S_->get_status()){
		do{ next_step(); }
		while(keepon());
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : MCSystem bad status (status="<<S_->get_status()<<")"<<std::endl; }
}

void MonteCarlo::run(unsigned int const& maxiter){
	time_.set();
	if(!S_->get_status()){
		Time chrono;
		unsigned int iter(0);
		unsigned int measures(0);
		do{
			S_->swap();
			ratio_ = my::norm_squared(S_->ratio());
			if( ratio_ > 1.0 || ratio_ > rnd_.get() ){
				S_->update();
				S_->measure_new_step();
				measures++;
			}
			S_->add_sample();
			if(iter%100000==0 && omp_get_thread_num()==0){ my::display_progress(time_.elapsed(),tmax_); }
		} while(keepon() && ++iter<maxiter);
#pragma omp critical(cout)
		std::cout<<"done "<<iter<<" steps in "<<chrono.elapsed()<<"s with "<<measures<<" measures"<<std::endl;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : MCSystem bad status (status="<<S_->get_status()<<")"<<std::endl; }
}
/*}*/

/*private methods*/
/*{*/
void MonteCarlo::next_step(){
	S_->swap();
	ratio_ = my::norm_squared(S_->ratio());
	if( ratio_ > 1.0 || ratio_ > rnd_.get() ){
		S_->update();
		S_->measure_new_step();
	}
	S_->add_sample();
}

bool MonteCarlo::keepon(){
	if(time_.limit_reached(tmax_)){ return false; }
	if(std::abs(S_->get_energy().get_x())>1e2){
		std::cerr<<__PRETTY_FUNCTION__<<" : simulation diverges (E="<<S_->get_energy().get_x()<<") => is restarted"<<std::endl;
		S_->reset_obs();
	}
	return true;
}
/*}*/
