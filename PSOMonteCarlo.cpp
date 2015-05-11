#include "PSOMonteCarlo.hpp"

PSOMonteCarlo::PSOMonteCarlo(Parseur& P):
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp")),
	tmax_(P.get<unsigned int>("tmax"))
{
	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);
	system_param_.set("wf",(in?in->read<std::string>():P.get<std::string>("wf")));
	system_param_.set("N",(in?in->read<unsigned int>():P.get<unsigned int>("N")));
	system_param_.set("m",(in?in->read<unsigned int>():P.get<unsigned int>("m")));
	system_param_.set("n",(in?in->read<unsigned int>():P.get<unsigned int>("n")));
	system_param_.set("bc",(in?in->read<int>():P.get<int>("bc")));
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
			sim->get_S()->set(); //to clear the binning
		} else {
			tmp_test = false;
			sim->create_S(&system_param_,&P->get_x());
		}
	}
	if(sim->is_created()){
		MonteCarlo mc(sim->get_S().get(),tmax_);
		mc.thermalize(tmp_test?10:1e6);
		mc.run();

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

void PSOMonteCarlo::refine(unsigned int const& Nrefine, double const& tol, unsigned int const& tmax){
	if(all_results_.size()){
		pso_info_.text("refine called with"+RST::nl_);
		pso_info_.item(my::tostring(Nrefine));
		pso_info_.item(my::tostring(tol));
		pso_info_.item(my::tostring(tmax)+RST::nl_);
		auto sort = [](MCSim const& a, MCSim const& b){ 
			return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
		};

		List<MCSim> best;
		while(all_results_.target_next()){ best.add_sort(all_results_.get_ptr(),sort); }
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
			MonteCarlo mc(sim->get_S().get(),tmax);
			while(sim->get_S()->get_energy().get_dx()>tol){
				mc.run();
				sim->get_S()->complete_analysis(1e-5);
			}
#pragma omp critical
			{
				std::cout<<i<<" sim refined "<<sim->get_S()->get_energy()<<std::endl;
			}
		}
	} else {
		std::cerr<<"void PSOMonteCarlo::refine(unsigned int const& Nrefine, double const& tol, unsigned int const& tmax) : there is no data"<<std::endl;
	}
	//best.set_free();
	//while(best.target_next()){
	//std::cout<<best.get().get_S()->get_energy()<<std::endl;
	//}
}

void PSOMonteCarlo::complete_analysis(double const& tol){
	all_results_.set_target();
	while ( all_results_.target_next() ){
		all_results_.get().get_S()->complete_analysis(tol);
	}
}

std::string PSOMonteCarlo::get_filename() const {
	std::string filename("PSO");
	filename += "-wf"+system_param_.get<std::string>("wf");
	filename += "-N" +my::tostring(system_param_.get<unsigned int>("N"));
	filename += "-m" +my::tostring(system_param_.get<unsigned int>("m"));
	filename += "-n" +my::tostring(system_param_.get<unsigned int>("n"));
	filename += "-bc"+my::tostring(system_param_.get<int>("bc"));
	filename += "_"+Time().date();
	return filename;
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
	gp+="plot '"+filename+".dat' u "+my::tostring(N+1)+":1:"+my::tostring(N+2)+" w xe notitle"+(N==1?"":",\\"); 
	for(unsigned int i(1);i<N;i++){
		gp+="     '"+filename+".dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe notitle"+(i==N-1?"":",\\");
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
		auto sort = [](MCSim const& a, MCSim const& b){ 
			return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
		};

		List<MCSim> best;
		while(all_results_.target_next()){ best.add_sort(all_results_.get_ptr(),sort); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(&system_param_); }
	} else {
		std::cerr<<"void PSOMonteCarlo::save(unsigned int const& nsave) : there is no data"<<std::endl;
	}
	//best.set_free();
	//while(best.target_next()){
	//std::cout<<best.get().get_S()->get_energy()<<std::endl;
	//}
}
