#include "VMCACiD.hpp"

VMCACiD::VMCACiD(VMCMinimization const& min, Vector<unsigned int> const& d):
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

		List<MCSim>::Node* target(get_best_target());

		Linux cmd;
		target->get()->display_results();
		cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
		std::cout<<target->get()->get_energy()<<std::endl;
		std::cout<<target->get()->get_param()<<std::endl;

		Vector<double> x(target->get()->get_param());
		sigma_.set(N_,0.05);
		Rand<double> rnd(-0.05,0.05);
		for(unsigned int i(0);i<m_->dof_;i++){
			idx_(i,0) = d(i);
			if(idx_(i,0)){
				idx_(i,1) = d(i)-1;
				xmean_(idx_(i,1)) = std::abs(x(i))+rnd();
			}
			idx_(i,2) = my::sign(x(i));
		}

		for(unsigned int i(0);i<m_->dof_;i++){
			x(i) = (idx_(i,0)?xmean_(idx_(i,1)):1.0)*(idx_(i,2)?idx_(i,2):1);
		}

		MCSim test(x);
		test.create_S(m_->s_);
		test.display_results();
		cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
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

	List<MCSim>::Node* target(get_best_target());

	Linux cmd;
	target->get()->display_results();
	cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
	std::cout<<target->get()->get_energy()<<std::endl;

	Vector<double> x(m_->dof_);
	for(unsigned int i(0);i<m_->dof_;i++){
		x(i) = (idx_(i,0)?xmean_(idx_(i,1)):1.0)*(idx_(i,2)?idx_(i,2):1);
	}
	MCSim test(x);
	test.create_S(m_->s_);
	test.display_results();

	cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
}

double VMCACiD::function(Vector<double> const& x){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0);i<m_->dof_;i++){
		param(i) = (idx_(i,0)?x(idx_(i,1)):1.0)*(idx_(i,2)?idx_(i,2):1.0);
	}
	Vector<unsigned int> obs(0,0);
	std::shared_ptr<MCSim> sim(evaluate_until_precision(param,obs,dEoE_,maxiter_));
	return sim->get_energy().get_x();
}

void VMCACiD::init(unsigned int const& maxiter, double const& dEoE, unsigned int const& tmax){
	m_->tmax_ = tmax;
	if(m_->tmax_ && idx_.ptr()){
		dEoE_ = dEoE;
		maxiter_ = maxiter;

		x_min_ = xmean_;
		bf_ = function(x_min_);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

void VMCACiD::run(unsigned int const& maxsteps){
	if(m_->tmax_ && idx_.ptr()){
		progress_ = 0;
		total_eval_ = 2*maxsteps;

		std::string msg("will perform a minimization with the ACiD method measuring " + my::tostring(total_eval_) + " new samples");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E/E="+my::tostring(dEoE_)) + " , maxiter="+my::tostring(maxiter_);
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		msg = "maximum time : "+my::tostring(total_eval_*m_->tmax_*maxiter_)+"s";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		ACiD::run(maxsteps);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

bool VMCACiD::keepon(bool const& improve_overall) const {
	if(!improve_overall){ std::cerr<<__PRETTY_FUNCTION__<<" : no imrovement detected, will continue anyway"<<std::endl; }
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
