#include "VMCMinimization.hpp"

/*{VMCMinimization*/
VMCMinimization::VMCMinimization(Parseur& P):
	time_(""),
	basename_(""),
	out_(NULL),
	m_(std::make_shared<Minimization>(P))
{
	set_time();
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

/*{public methods*/
void VMCMinimization::refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax){
	if(m_->samples_list_.size()){
		std::cout<<"#######################"<<std::endl;
		std::string msg("refine ("+my::tostring(Nrefine)+","+my::tostring(convergence_criterion)+","+my::tostring(tmax)+")");
		std::cout<<"#"<<msg<<std::endl;
		m_->pso_info_.item(msg);

		List<MCSim> best;
		while(m_->samples_list_.target_next()){ best.add_sort(m_->samples_list_.get_ptr(),MCSim::compare); }
		unsigned int N(m_->samples_list_.size());
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
	m_->samples_list_.set_target();
	while ( m_->samples_list_.target_next() ){
		m_->samples_list_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	IOFiles out(get_filename()+".jdbin",true);
	m_->save(out);
}

void VMCMinimization::save_best(unsigned int const& nsave){
	if(m_->samples_list_.size()){
		List<MCSim> best;
		while(m_->samples_list_.target_next()){ best.add_sort(m_->samples_list_.get_ptr(),MCSim::compare); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(&m_->system_param_); }
	} else {
		std::cerr<<"void VMCMinimization::save(unsigned int const& nsave) : there is no data"<<std::endl;
	}
}

void VMCMinimization::print() const {
	std::cout<<"Print whole history ("<< m_->samples_list_.size()<<")"<<std::endl;
	while( m_->samples_list_.target_next() ){
		std::cout<<m_->samples_list_.get_ptr()<<" ";
		m_->samples_list_.get().print();
		std::cout<<std::endl;
	}
}

void VMCMinimization::plot() const {
	std::string filename("PS"+get_filename());
	IOFiles data(filename+".dat",true);
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		data<<m_->samples_list_.get().get_param()<<" "<<m_->samples_list_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	m_->samples_list_.target_next();
	Gnuplot gp("./",filename);
	unsigned int N(m_->samples_list_.get().get_param().size());
	for(unsigned int i(0);i<N;i++){
		gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe notitle"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true);
}
/*}*/

/*{protected methods*/
std::shared_ptr<MCSim> VMCMinimization::evaluate(Vector<double> const& param){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	bool tmp_test;
#pragma omp critical(samples_list_)
	{
		if(m_->samples_list_.find_sorted(sim,MCSim::cmp_for_merge)){ 
			tmp_test = true;
			sim->copy_S(m_->samples_list_.get().get_S()); 
		} else {
			tmp_test = false;
			sim->create_S(&m_->system_param_);
		}
		m_->samples_list_.set_target();
	}
	if(sim->is_created()){
		sim->run(tmp_test?10:1e6,m_->tmax_);
#pragma omp critical(samples_list_)
		{
			if(m_->samples_list_.find_sorted(sim,MCSim::cmp_for_merge)){ 
				m_->samples_list_.merge_with_target(sim,MCSim::merge); 
				sim = m_->samples_list_.get_ptr();
			} else {
				m_->samples_list_.add_after_target(sim); 
			}
			m_->samples_list_.set_target();
		}
		sim->free_memory();
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
	tmax_(P.get<unsigned int>("tmax"))
{
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#creating Minimization"<<std::endl;

	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);

	wf_      =(in?in->read<std::string>() :P.get<std::string>("wf"));
	N_       =(in?in->read<unsigned int>():P.get<unsigned int>("N"));
	m_       =(in?in->read<unsigned int>():P.get<unsigned int>("m"));
	n_       =(in?in->read<unsigned int>():P.get<unsigned int>("n"));
	bc_      =(in?in->read<int>()         :P.get<int>("bc"));
	Nfreedom_=(in?in->read<unsigned int>():P.get<unsigned int>("Nfreedom"));
	ps_ = new Vector<double>[Nfreedom_];
	for(unsigned int i(0);i<Nfreedom_;i++){
		ps_[i] = (in?in->read<Vector<double> >():P.get<std::vector<double> >("ps"+my::tostring(i))); 
	}
	if(in){
		std::string msg1("loading samples from "+in->get_filename());
		std::cout<<"#"+msg1<<std::flush;
		Time chrono;
		unsigned int size(in->read<int>());
		while(size--){ samples_list_.add_end(std::make_shared<MCSim>(*in)); }

		std::string msg2(" ("+my::tostring(samples_list_.size())+" samples loaded in "+my::tostring(chrono.elapsed())+"s)");
		std::string header(in->get_header());
		header.erase(0,header.find(".. end_of_saved_variables\n",0)+28);
		pso_info_.text(header);
		pso_info_.title("Minimization",'>');
		pso_info_.item(msg1+msg2);
		std::cout<<msg2<<std::endl;

		delete in;
		in = NULL;
	} else {
		std::string msg("no samples loaded");
		std::cout<<"#"+msg<<std::endl;
		pso_info_.title("Minimization",'>');
		pso_info_.item(msg);
	}

	system_param_.set("wf",wf_);
	system_param_.set("N",N_);
	system_param_.set("m",m_);
	system_param_.set("n",n_);
	system_param_.set("bc",bc_);
}

VMCMinimization::Minimization::~Minimization(){
	std::cerr<<pso_info_.get()<<std::endl;
	if(ps_){ delete[] ps_; }
}

bool VMCMinimization::Minimization::within_limit(Vector<double> const& x){
	for(unsigned int i(0);i<Nfreedom_;i++){
		if(x(i)<ps_[i](0) || x(i)>ps_[i].back()) { return false; }
	}
	return true;
}

void VMCMinimization::Minimization::save(IOFiles& out) const {
	out.write("wf",wf_);
	out.write("N", N_);
	out.write("m", m_);
	out.write("n", n_);
	out.write("bc",bc_);
	out.write("Nfreedom",Nfreedom_);
	for(unsigned int i(0);i<Nfreedom_;i++){ out<<ps_[i]; }
	out.write("# samples",samples_list_.size());
	out.add_header()->nl();
	out.add_header()->comment("end_of_saved_variables");
	out.add_header()->text(pso_info_.get());
	while(samples_list_.target_next()){ samples_list_.get().write(out); }
}
/*}*/
