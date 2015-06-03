#include "VMCMinimization.hpp"

/*{VMCMinimization*/
VMCMinimization::VMCMinimization(Parseur& P):
	time_(""),
	basename_(""),
	out_(NULL),
	m_(std::make_shared<Minimization>(P))
{
	basename_ += "-" + m_->wf_;
	basename_ += "-N"  + my::tostring(m_->N_);
	basename_ += "-m"  + my::tostring(m_->m_);
	basename_ += "-n"  + my::tostring(m_->n_);
	basename_ += "-bc" + my::tostring(m_->bc_);
}

VMCMinimization::VMCMinimization(VMCMinimization const& vmcm, std::string const& prefix):
	time_(""),
	basename_(prefix+vmcm.basename_),
	out_(NULL),
	m_(vmcm.m_)
{}

/*{Public methods*/
void VMCMinimization::refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax){
	if(m_->all_results_.size()){
		std::cout<<"#######################"<<std::endl;
		std::cout<<"#refine called with Nrefine="<<Nrefine<<" convergence_criterion="<<convergence_criterion<<" tmax="<<tmax<<std::endl;
		m_->pso_info_.text("refine called with"+RST::nl_);
		m_->pso_info_.item(my::tostring(Nrefine));
		m_->pso_info_.item(my::tostring(convergence_criterion));
		m_->pso_info_.item(my::tostring(tmax)+RST::nl_);

		List<MCSim> best;
		while(m_->all_results_.target_next()){ best.add_sort(m_->all_results_.get_ptr(),MCSim::compare); }
		unsigned int N(m_->all_results_.size());
		N  = (N<Nrefine?N:Nrefine);
		best.set_target();
#pragma omp parallel for schedule(dynamic,1)
		for(unsigned int i=0;i<N;i++){
			std::shared_ptr<MCSim> sim;
#pragma omp critical
			{
				best.target_next();
				sim = best.get_ptr();
			}
			while(!sim->check_conv(convergence_criterion)) { sim->run(0,tmax); }
			sim->complete_analysis(convergence_criterion);

#pragma omp critical
			{
				std::cout<<i<<" sim refined "<<sim->get_S()->get_energy()<<std::endl;
			}
		}
	} else {
		std::cerr<<"void PSO::refine(unsigned int const& Nrefine, double const& converged_criterion, unsigned int const& tmax) : there is no data"<<std::endl;
	}
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#complete_analysis called with convergence_criterion="<<convergence_criterion<<std::endl;
	m_->all_results_.set_target();
	while ( m_->all_results_.target_next() ){
		m_->all_results_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	IOFiles out(get_filename()+".jdbin",true);
	out.write("wf",m_->system_param_.get<std::string>("wf"));
	out.write("N", m_->system_param_.get<unsigned int>("N"));
	out.write("m", m_->system_param_.get<unsigned int>("m"));
	out.write("n", m_->system_param_.get<unsigned int>("n"));
	out.write("bc",m_->system_param_.get<int>("bc"));
	out.write("#", m_->all_results_.size());
	out.add_header()->nl();
	out.add_header()->text(m_->pso_info_.get());
	while(m_->all_results_.target_next()){ m_->all_results_.get().write(out); }
}

void VMCMinimization::save_best(unsigned int const& nsave){
	if(m_->all_results_.size()){
		List<MCSim> best;
		while(m_->all_results_.target_next()){ best.add_sort(m_->all_results_.get_ptr(),MCSim::compare); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(&m_->system_param_); }
	} else {
		std::cerr<<"void VMCMinimization::save(unsigned int const& nsave) : there is no data"<<std::endl;
	}
}

/*{Virtual methods*/
void VMCMinimization::set_ps(unsigned int const& i, Vector<double> const& ps){
	if(i<m_->Nfreedom_){
		if(!m_->ps_){ m_->ps_ = new Vector<double>[m_->Nfreedom_]; }
		m_->ps_[i] = ps;
	} else {
		std::cerr<<"void Minimization::set_x(unsigned int const& i, Vector<double> const& x) : i>=m_->Nfreedom"<<std::endl;
	}
}

void VMCMinimization::print() const {
	std::cout<<"Print whole history ("<< m_->all_results_.size()<<")"<<std::endl;
	while( m_->all_results_.target_next() ){
		std::cout<<m_->all_results_.get_ptr()<<" ";
		m_->all_results_.get().print();
		std::cout<<std::endl;
	}
}
/*}*/
/*}*/

/*{Protected methods*/
std::shared_ptr<MCSim> VMCMinimization::evaluate(Vector<double> const& param){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	bool tmp_test;
#pragma omp critical(all_results_)
	{
		if(m_->all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
			tmp_test = true;
			sim->copy_S(m_->all_results_.get().get_S()); 
		} else {
			tmp_test = false;
			sim->create_S(&m_->system_param_);
		}
	}
	if(sim->is_created()){
		sim->run(tmp_test?10:1e6,m_->tmax_);
#pragma omp critical(all_results_)
		{
			if(m_->all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
				m_->all_results_.fuse_with_target(sim,MCSim::fuse); 
				sim = m_->all_results_.get_ptr();
			} else {
				m_->all_results_.add_after_target(sim); 
			}
		}
		return sim;
	} else {
		std::cerr<<"bool Minimization::evaluate(Vector<double> const& param) : not valid parameter : "<<param<<std::endl;
		return NULL;
	}
}
/*}*/
/*}*/

/*{Minimization*/
VMCMinimization::Minimization::Minimization(Parseur& P):
	Nfreedom_(P.get<unsigned int>("Nfreedom")),
	tmax_(P.get<unsigned int>("tmax")),
	ps_(NULL)
{
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#creating Minimization"<<std::endl;
	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);

	wf_=(in?in->read<std::string>():P.get<std::string>("wf"));
	N_ =(in?in->read<unsigned int>():P.get<unsigned int>("N"));
	m_ =(in?in->read<unsigned int>():P.get<unsigned int>("m"));
	n_ =(in?in->read<unsigned int>():P.get<unsigned int>("n"));
	bc_=(in?in->read<int>():P.get<int>("bc"));

	system_param_.set("wf",wf_);
	system_param_.set("N",N_);
	system_param_.set("m",m_);
	system_param_.set("n",n_);
	system_param_.set("bc",bc_);

	if(in){
		std::cout<<"#loading";
		Time chrono;
		unsigned int size(in->read<int>());
		while(size--){ all_results_.add_end(std::make_shared<MCSim>(*in)); }
		pso_info_.text("loads data from "+in->get_filename()+RST::nl_);

		delete in;
		in = NULL;
		std::cout<<" ("<<all_results_.size()<<" samples loaded "<<chrono.elapsed()<<"s)"<<std::endl;
	} else {
		std::cout<<"#no samples loaded"<<std::endl;
	}
}

VMCMinimization::Minimization::~Minimization(){
	if(ps_){ delete[] ps_; }
}

bool VMCMinimization::Minimization::within_limit(Vector<double> const& x){
	for(unsigned int i(0);i<Nfreedom_;i++){
		if(x(i)<ps_[i](0) || x(i)>ps_[i].back()) { return false; }
	}
	return true;
}
/*}*/
