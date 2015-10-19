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
			ratio_ = S_->ratio(true);
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
	}
}
/*}*/

/*private methods*/
/*{*/
void MonteCarlo::next_step(){
	S_->swap();
	ratio_ = S_->ratio(true);
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
		S_->clear_measurments();
	}
	return true;
}
/*}*/
