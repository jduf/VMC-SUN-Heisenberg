#include "VMCMinimization.hpp"

VMCMinimization::VMCMinimization(Parseur& P, std::string const& prefix):
	Nfreedom_(P.get<unsigned int>("Nfreedom")),
	tmax_(P.get<unsigned int>("tmax")),
	ps_(NULL),
	basename_(prefix)
{
	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);
	std::string wf(in?in->read<std::string>():P.get<std::string>("wf"));
	unsigned int N(in?in->read<unsigned int>():P.get<unsigned int>("N"));
	unsigned int m(in?in->read<unsigned int>():P.get<unsigned int>("m"));
	unsigned int n(in?in->read<unsigned int>():P.get<unsigned int>("n"));
	int bc(in?in->read<int>():P.get<int>("bc"));

	system_param_.set("wf",wf);
	system_param_.set("N",N);
	system_param_.set("m",m);
	system_param_.set("n",n);
	system_param_.set("bc",bc);

	basename_ += "-wf"+system_param_.get<std::string>("wf");
	basename_ += "-N" +my::tostring(N);
	basename_ += "-m" +my::tostring(m);
	basename_ += "-n" +my::tostring(n);
	basename_ += "-bc"+my::tostring(bc);

	if(in){
		unsigned int size(in->read<int>());
		while(size--){ all_results_.add_end(std::make_shared<MCSim>(*in)); }
		pso_info_.text("loads data from "+in->get_filename()+RST::nl_);

		delete in;
		in = NULL;
	}
}

VMCMinimization::~VMCMinimization(){
	if(ps_){ delete[] ps_; }
}

void VMCMinimization::set_ps(unsigned int const& i, Vector<double> const& ps){
	if(i<Nfreedom_){
		if(!ps_){ ps_ = new Vector<double>[Nfreedom_]; }
		ps_[i] = ps;
	} else {
		std::cerr<<"void VMCMinimization::set_x(unsigned int const& i, Vector<double> const& x) : i>=Nfreedom"<<std::endl;
	}
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#complete_analysis called with convergence_criterion="<<convergence_criterion<<std::endl;
	all_results_.set_target();
	while ( all_results_.target_next() ){
		all_results_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	IOFiles out(get_filename()+".jdbin",true);
	out.write("wf",system_param_.get<std::string>("wf"));
	out.write("N", system_param_.get<unsigned int>("N"));
	out.write("m", system_param_.get<unsigned int>("m"));
	out.write("n", system_param_.get<unsigned int>("n"));
	out.write("bc",system_param_.get<int>("bc"));
	out.write("#", all_results_.size());
	out.add_header()->nl();
	out.add_header()->text(pso_info_.get());
	while(all_results_.target_next()){ all_results_.get().write(out); }
}

void VMCMinimization::refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax){
	if(all_results_.size()){
		std::cout<<"#######################"<<std::endl;
		std::cout<<"#refine called with Nrefine="<<Nrefine<<" convergence_criterion="<<convergence_criterion<<" tmax="<<tmax<<std::endl;
		pso_info_.text("refine called with"+RST::nl_);
		pso_info_.item(my::tostring(Nrefine));
		pso_info_.item(my::tostring(convergence_criterion));
		pso_info_.item(my::tostring(tmax)+RST::nl_);

		List<MCSim> best;
		while(all_results_.target_next()){ best.add_sort(all_results_.get_ptr(),VMCMinimization::sort_per_energy); }
		unsigned int N(all_results_.size());
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
		std::cerr<<"void VMCPSO::refine(unsigned int const& Nrefine, double const& converged_criterion, unsigned int const& tmax) : there is no data"<<std::endl;
	}
}

void VMCMinimization::move(VMCMinimization* min){
	all_results_.move(min->all_results_);
	pso_info_.text(min->pso_info_.get());
	min->pso_info_.text("");
}

std::shared_ptr<MCSim> VMCMinimization::compute_vmc(Vector<double> const& param){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	bool tmp_test;
#pragma omp critical(all_results_)
	{
		if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
			tmp_test = true;
			sim->copy_S(all_results_.get().get_S()); 
		} else {
			tmp_test = false;
			sim->create_S(&system_param_);
		}
	}
	if(sim->is_created()){
		sim->run(tmp_test?10:1e6,tmax_);
#pragma omp critical(all_results_)
		{
			if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
				all_results_.fuse_with_target(sim,MCSim::fuse); 
				sim = all_results_.get_ptr();
			} else {
				all_results_.add_after_target(sim); 
			}
		}
		return sim;
	} else {
		std::cerr<<"bool VMCMinimization::compute_vmc(Vector<double> const& param) : not valid parameter : "<<param<<std::endl;
		return NULL;
	}
}
