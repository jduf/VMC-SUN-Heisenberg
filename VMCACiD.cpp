#include "VMCACiD.hpp"

VMCACiD::VMCACiD(VMCMinimization const& min, Vector<unsigned int> const& d, Vector<double> const& param):
	VMCMinimization(min,"ACiD"),
	ACiD(d.max(),5.0),
	idx_(m_->dof_,3,0)
{
	if(m_->dof_ != d.size() && d.max() > m_->dof_){
		std::cerr<<__PRETTY_FUNCTION__<<" : wrong number of free parameters"<<std::endl;
		std::cerr<<__PRETTY_FUNCTION__<<"m_->dof_="<<m_->dof_<<"!="<<d.size()<<"=d.size() || m_->dof_="<<m_->dof_<<">"<<d.max()<<"=d.max()"<<std::endl;
		idx_.set();
	} else {
		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCACiD");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		Vector<double> x;
		if(param.size()){ x = param; }
		else {
			List<MCSim>::Node* target(get_best_target());
			x = target->get()->get_param();
		}

		sigma_.set(N_,0.05);
		RandGaussian rnd(0,0.025);
		for(unsigned int i(0);i<m_->dof_;i++){
			idx_(i,0) = d(i);
			if(idx_(i,0)){
				idx_(i,1) = d(i)-1;
				xmean_(idx_(i,1)) = std::abs(x(i))+rnd();
			}
			idx_(i,2) = my::sign(x(i));
		}
	}
}

VMCACiD::VMCACiD(VMCMinimization const& min, IOFiles& in):
	VMCMinimization(min,"ACiD"),
	ACiD(in),
	idx_(in)
{
	std::cout<<"#######################"<<std::endl;
	std::string msg("VMCACiD");
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.title(msg,'-');
}

double VMCACiD::function(Vector<double> const& x){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0);i<m_->dof_;i++){
		param(i) = (idx_(i,0)?x(idx_(i,1)):1.0)*(idx_(i,2)?idx_(i,2):1.0);
	}
	std::shared_ptr<MCSim> sim(evaluate_until_precision(param,0,dEoE_,maxiter_));
	return sim->get_energy().get_x();
}

void VMCACiD::run(double const& dEoE, unsigned int const& maxiter, unsigned int const& tmax, unsigned int const& maxsteps, std::string const& save_in){
	m_->tmax_ = tmax;
	dEoE_ = dEoE;
	maxiter_ = maxiter;

	if(m_->tmax_ && idx_.ptr()){
		std::string msg(RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E/E="+my::tostring(dEoE_)) + ", maxiter="+my::tostring(maxiter_));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(!x_min_.ptr()){
			msg = "initialise the starting point";
			std::cout<<"#"<<msg<<std::endl;
			m_->info_.item(msg);

			progress_ = 0;
			total_eval_ = 1;
			x_min_ = xmean_;
			bf_ = function(x_min_);
		}

		progress_ = 0;
		total_eval_ = 2*maxsteps;

		msg = "will perform a minimization with the ACiD method measuring " + my::tostring(total_eval_) + " new samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		msg = "maximum time : "+my::tostring(total_eval_*m_->tmax_*maxiter_)+"s";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		ACiD::run(maxsteps);
		save(save_in);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

bool VMCACiD::keepon(bool const& improve_overall) const {
	if(!improve_overall){ std::cerr<<__PRETTY_FUNCTION__<<" : no improvement detected, will continue anyway"<<std::endl; }
	return true;
}

void VMCACiD::save(std::string dirname) const {
	my::ensure_trailing_slash(dirname);
	dirname += get_path();
	set_time();
	Linux().mkpath(dirname.c_str());
	IOFiles out(dirname+get_filename()+".jdbin",true,false);
	ACiD::save(out);
	out<<idx_;
}

void VMCACiD::display_param_and_xmean(Vector<double> const& param) const {
	List<MCSim>::Node* target(NULL);
	Vector<double> x(m_->dof_);
	if(param.size()){
		std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
		if(m_->samples_.find_in_sorted_list(sim,target,MCSim::sort_for_merge)){
			x = param;
		} else {
			std::cerr<<__PRETTY_FUNCTION__<<" : can't find sample with parameter "<<param<<std::endl;
		}
	} else {
		target = get_best_target();
		x = target->get()->get_param();
	}

	Linux cmd;
	if(target){
		target->get()->display_results();
		cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
		std::cout<<target->get()->get_energy()<<std::endl;
		std::cout<<target->get()->get_param()<<std::endl;
	}

	for(unsigned int i(0);i<m_->dof_;i++){
		x(i) = (idx_(i,0)?xmean_(idx_(i,1)):1.0)*(idx_(i,2)?idx_(i,2):1);
	}

	MCSim test(x);
	test.create_S(m_->s_);
	test.display_results();
	cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
}

void VMCACiD::test(){
	Vector<double> x(m_->dof_);
	for(unsigned int i(0);i<m_->dof_;i++){
		x(i) = (idx_(i,0)?xmean_(idx_(i,1)):1.0)*(idx_(i,2)?idx_(i,2):1);
	}

	std::cout<<function(x)<<std::endl;
}
