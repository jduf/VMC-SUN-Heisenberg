#include "VMCSystematic.hpp"

VMCSystematic::VMCSystematic(VMCMinimization const& m):
	VMCMinimization(m,"SYS")
{}

/*{public methods*/
void VMCSystematic::run(int const& nobs, double const& dE, unsigned int const& maxiter){
	nobs_ = nobs;
	dE_ = dE;
	maxiter_ = maxiter;
	Vector<unsigned int> idx(m_->dof_,0);
	while(go_through_parameter_space(m_->ps_,idx,0,0,&VMCSystematic::evaluate));
}

void VMCSystematic::plot(){
	IOFiles data(get_path()+get_filename()+".dat",true);
	m_->samples_.set_target();
	Vector<double> param;
	while(m_->samples_.target_next()){
		param = m_->samples_.get().get_param(); 
		data<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_MCS()->get_energy()<<IOFiles::endl;
	}

	if(m_->dof_==1){
		Gnuplot gp(get_path(),get_filename());
		gp+="plot '"+get_filename()+".dat' u 1:2:3 w e notitle";
		gp.save_file();
		gp.create_image(true,true);
	} else {
		Gnuplot gp(get_path(),get_filename());
		gp+="plot '"+get_filename()+".dat' u  w e notitle";
		gp.save_file();
	}
}
/*}*/

bool VMCSystematic::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSystematic::*f)(Vector<double>*, Vector<unsigned int> const&)){
	(this->*f)(x,idx);

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<m_->dof_;i++){
			if(idx(i) == x[i].size()){
				if(i+1 == m_->dof_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	} else {
		if(idx(0) == max0){
			if(1 == m_->dof_){ return false; }
			idx(0) = min0;
			idx(1)++;
		}
		for(unsigned int i(1);i<m_->dof_;i++){
			if(idx(i) == x[i].size()){
				if(i+1 == m_->dof_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return go_through_parameter_space(x,idx,min0,max0,f);
}

void VMCSystematic::evaluate(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0); i<m_->dof_;i++){ param(i) = x[i](idx(i)); }
	evaluate_until_precision(param,nobs_,dE_,maxiter_);
}
