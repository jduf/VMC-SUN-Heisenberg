#include "VMCACiD.hpp"

VMCACiD::VMCACiD(VMCMinimization const& min, Vector<unsigned int> const& d):
	VMCMinimization(min,"ACiD"),
	ACiD(d.max(),Vector<double>(m_->dof_,-0.05),Vector<double>(m_->dof_,0.05),5.0),
	idx_(m_->dof_,3,0)
{
	if(m_->dof_ != d.size()){
		std::cerr<<__PRETTY_FUNCTION__<<" : wrong number of element to set for d"<<std::endl;
	} else {
		List<MCSim>::Node* target(get_best_target());

		Linux cmd;
		target->get()->display_results();
		cmd("~/divers/linux/scripts/display-image-terminal.bash "+cmd.pwd()+"tmp.png min",false);
		std::cout<<target->get()->get_energy()<<std::endl;


		Vector<double> x(target->get()->get_param());
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

void VMCACiD::run(unsigned int const& maxsteps, unsigned int const& maxiter, double const& dEoE){
	if(m_->tmax_){
		progress_ = 0;
		total_eval_ = 2*maxsteps+1;
		dEoE_ = dEoE;
		maxiter_ = maxiter;

		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCACiD");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		msg = "will perform a minimization with the ACiD method measuring " + my::tostring(total_eval_) + " new samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E/E="+my::tostring(dEoE)) + " , maxiter="+my::tostring(maxiter_);
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		msg = "maximum time : "+my::tostring(total_eval_*m_->tmax_*maxiter)+"s";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		ACiD::run(maxsteps);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

bool VMCACiD::stop(bool const& improve_overall) const {
	if(!improve_overall){ std::cerr<<__PRETTY_FUNCTION__<<" : no imrovement detected, will continue anyway"<<std::endl; }
	return false;
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
