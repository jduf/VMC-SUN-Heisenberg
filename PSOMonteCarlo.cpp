#include "PSOMonteCarlo.hpp"

PSOMonteCarlo::PSOMonteCarlo(Parseur& P):
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp")),
	tmax_(P.get<unsigned int>("tmax")),
	basename_("PSO")
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

void PSOMonteCarlo::init(bool const& clear_particle_history, bool const& create_particle_history){
	pso_info_.title("New run",'-');
	if(clear_particle_history){ 
		pso_info_.text("clear history"+RST::nl_);
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
		}
	}
	time_ = Time().date();

	track_particles_ = new IOFiles(get_filename()+"-ERR.dat",true);
	init_PSO(100); 
	if(clear_particle_history && create_particle_history && all_results_.size()){
		pso_info_.text("set history"+RST::nl_);
		unsigned int size(0);
		while(all_results_.target_next()){
			if(Optimization::within_limit(all_results_.get().get_param())){ size++; }
		}
		unsigned int Npp(size/Nparticles_);
		unsigned int s;
		std::shared_ptr<MCParticle> P;
		for(unsigned int p(0);p<Nparticles_;p++){
			if(p==Nparticles_-size%Nparticles_){ Npp++; }
			P = std::dynamic_pointer_cast<MCParticle>(particle_[p]);
			s = 0;
			while( s!=Npp && all_results_.target_next()){
				if(Optimization::within_limit(all_results_.get().get_param())){
					s++;
					P->add_to_history(all_results_.get_ptr());
				}
			}
			P->select_new_best();
		}
	}
}

bool PSOMonteCarlo::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> P(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(P->get_x()));
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
				tmp_test = P->update(all_results_.get_ptr()); 
			} else {
				all_results_.add_after_target(sim); 
				tmp_test = P->update(sim); 
			}
		}
		return tmp_test;
	} else {
		std::cerr<<"bool PSOMonteCarlo::evaluate(unsigned int const& p) : not valid parameter : "<<particle_[p]->get_x()<<std::endl;
		return false;
	}
}

void PSOMonteCarlo::refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax){
	if(all_results_.size()){
		pso_info_.text("refine called with"+RST::nl_);
		pso_info_.item(my::tostring(Nrefine));
		pso_info_.item(my::tostring(convergence_criterion));
		pso_info_.item(my::tostring(tmax)+RST::nl_);

		List<MCSim> best;
		while(all_results_.target_next()){ best.add_sort(all_results_.get_ptr(),PSOMonteCarlo::sort_per_energy); }
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
		std::cerr<<"void PSOMonteCarlo::refine(unsigned int const& Nrefine, double const& converged_criterion, unsigned int const& tmax) : there is no data"<<std::endl;
	}
}

void PSOMonteCarlo::complete_analysis(double const& converged_criterion){
	all_results_.set_target();
	while ( all_results_.target_next() ){
		all_results_.get().complete_analysis(converged_criterion);
	}

	if(track_particles_){
		Gnuplot gp("./",get_filename()+"-ERR");
		gp+="ERRFILE='"+get_filename()+"-ERR.dat'";
		Vector<double> m(get_bx());
		std::string space(16+my::tostring(Nparticles_).size(),' ');
		for(unsigned int i(0);i<Optimization::get_Nfreedom();i++){ gp+="a"+my::tostring(i)+"="+my::tostring(m(i)); }
		gp+="";
		gp+="set multiplot layout 2,1";
		for(unsigned int i(0);i<Optimization::get_Nfreedom();i++){
			gp+=std::string(!i?"plot":"    ")+" for [p=0:"+my::tostring(Nparticles_-1)+"] ERRFILE u 1:($2==p? ($"+my::tostring(4+i)+"-a"+my::tostring(i)+"):1/0) pt 13 lc "+my::tostring(i)+" notitle,\\";
		}
		for(unsigned int i(0);i<Optimization::get_Nfreedom();i++){
			gp+=space+"ERRFILE u 1:($2==-1? ($"+my::tostring(4+i)+"-a"+my::tostring(i)+"):1/0) pt 4 lc "+my::tostring(i)+" t sprintf('%f',a"+my::tostring(i)+")"+(i==Optimization::get_Nfreedom()-1?"":",\\");
		}
		gp+="";
		gp+="plot for [p=0:"+my::tostring(Nparticles_-1)+"] ERRFILE u 1:($2==p? $3:1/0) pt 13 lc 1  notitle,\\";
		gp+=space+"ERRFILE u 1:($2==-1?$3:1/0) pt 4 lc 1 notitle";
		gp+="unset multiplot";
		gp.save_file();

		delete track_particles_;
		track_particles_ = NULL;
	}
}

void PSOMonteCarlo::plot() const {
	std::string filename(get_filename());
	IOFiles data(filename+".dat",true);
	all_results_.set_target();
	while(all_results_.target_next()){
		data<<all_results_.get().get_param()<<all_results_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	all_results_.target_next();
	Gnuplot gp("./",filename);
	unsigned int N(all_results_.get().get_param().size());
	for(unsigned int i(0);i<N;i++){
		gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe notitle"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true);
}

void PSOMonteCarlo::print() const {
	Swarm::print();
	std::cout<<"Print whole history ("<< all_results_.size()<<")"<<std::endl;
	while( all_results_.target_next() ){
		std::cout<<all_results_.get_ptr()<<" ";
		all_results_.get().print();
		std::cout<<std::endl;
	}
}

void PSOMonteCarlo::save() const {
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

void PSOMonteCarlo::save(unsigned int const& nsave){
	if(all_results_.size()){
		List<MCSim> best;
		while(all_results_.target_next()){ best.add_sort(all_results_.get_ptr(),PSOMonteCarlo::sort_per_energy); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(&system_param_); }
	} else {
		std::cerr<<"void PSOMonteCarlo::save(unsigned int const& nsave) : there is no data"<<std::endl;
	}
}
