#include "VMCMinimization.hpp"

/*{VMCMinimization*/
VMCMinimization::VMCMinimization(Parseur& P):
	time_(""),
	basename_(""),
	prefix_("MIN"),
	out_(NULL),
	m_(std::make_shared<Minimization>())
{
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#creating VMCMinimization"<<std::endl;
	basename_ = m_->set(P);
	set_time();
	if(m_->s_->get_status() != 3 || P.locked()){
		std::cout<<m_->s_->get_status()<<std::endl;
		m_.reset();
		std::cerr<<"VMCMinimization::VMCMinimization(Parseur& P) : something went wrong"<<std::endl;
	}
}

VMCMinimization::VMCMinimization(VMCMinimization const& vmcm, std::string const& prefix):
	time_(""),
	basename_(vmcm.basename_),
	prefix_(prefix),
	out_(NULL),
	m_(vmcm.m_)
{}

/*{public methods*/
void VMCMinimization::refine(){
	double E;
	double dE(0.1);
	while(dE>1e-5){
		E=0;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_S()->get_energy().get_x()<E){
				E = m_->samples_list_.get().get_S()->get_energy().get_x();
			}
		}
		E += 1.5*dE;
		refine(E,dE);
		dE /= 2;
		if(dE<1e-3){ m_->tmax_ *= 2; }
	}
}

void VMCMinimization::refine(double const& E, double const& dE){
	if(m_->samples_list_.size()){
		List<MCSim> best;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_S()->get_energy().get_x()<E){
				best.add_end(m_->samples_list_.get_ptr()); 
			}
		}
		unsigned int N(best.size());

		std::cout<<"#######################"<<std::endl;
		std::string msg("refine "+my::tostring(N)+" samples ("+my::tostring(E)+","+my::tostring(dE)+","+my::tostring(m_->tmax_)+")");
		std::cout<<"#"<<msg<<std::endl;
		m_->pso_info_.item(msg);

		/*because there is no default MCSim constructor*/
		Vector<double> tmp(m_->Nfreedom_);
#pragma omp parallel for schedule(dynamic,10)
		for(unsigned int i=0;i<N;i++){
			MCSim sim(tmp);
#pragma omp critical
			{
				best.target_next();
				/*need to do call copy_S because best.get().get_S() has an
				 * undefined Ainv_, EVec_[i>0] ... due to the free_memory()
				 * call*/
				sim.copy_S(best.get().get_S());
			}
			while(!sim.check_conv(1e-5) || sim.get_S()->get_energy().get_dx()>dE) { sim.run(0,m_->tmax_); }
			sim.complete_analysis(dE);
		}
	} else {
		std::cerr<<"void VMCMinimization::refine(unsigned int const& Nrefine, double const& converged_criterion, unsigned int const& tmax) : there is no data"<<std::endl;
	}
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#complete_analysis called with convergence_criterion="<<convergence_criterion<<std::endl;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		m_->samples_list_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	IOFiles out(get_filename()+".jdbin",true);
	m_->save(out);
	out<<basename_;
}

void VMCMinimization::save_best(unsigned int const& nsave){
	if(m_->samples_list_.size()){
		List<MCSim> best;
		while(m_->samples_list_.target_next()){ best.add_sort(m_->samples_list_.get_ptr(),MCSim::compare); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(m_->s_); }
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
	std::string filename(get_filename());
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
			sim->create_S(m_->s_);
		}
		m_->samples_list_.set_target();
	}
	if(sim->is_created()){
		sim->set_observable(0);
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
std::string VMCMinimization::Minimization::set(Parseur& P){
	tmax_ = P.get<unsigned int>("tmax");

	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);

	Vector<unsigned int> ref(in?in->read<Vector<unsigned int> >():CreateSystem::get_ref(P.get<std::string>("wf")));
	unsigned int         N  (in?in->read<unsigned int>()         :P.get<unsigned int>("N"));
	unsigned int         m  (in?in->read<unsigned int>()         :P.get<unsigned int>("m"));
	unsigned int         n  (in?in->read<unsigned int>()         :P.get<unsigned int>("n"));
	unsigned int         bc (in?in->read<int>()                  :P.get<int>("bc"));
	Vector<unsigned int> M  (in?in->read<Vector<unsigned int> >():P.get<std::vector<unsigned int>>("M"));
	Vector<double>       J  (in?in->read<Vector<double> >()      :P.get<std::vector<double> >("Jp"));
	Nfreedom_             = (in?in->read<unsigned int>()         :P.get<unsigned int>("Nfreedom"));
	s_ = new System(ref,N,m,n,bc,M,J);
	ps_= new Vector<double>[Nfreedom_];

	/*!the next block is required to compute the J setup*/
	Vector<double> tmp(Nfreedom_,1);
	CreateSystem cs(s_);
	cs.set_param(NULL,&tmp);
	cs.init();
	cs.set_bonds(s_);

	std::string basename;
	ps_size_ = 1;
	if(in){
		for(unsigned int i(0);i<Nfreedom_;i++){
			ps_[i] = in->read<Vector<double> >(); 
			ps_size_ *= ps_[i].size();
		}

		std::string msg("loading samples from "+in->get_filename());
		std::cout<<"#"+msg<<std::flush;
		Time chrono;
		unsigned int size(in->read<int>());
		while(size--){ samples_list_.add_end(std::make_shared<MCSim>(*in)); }
		(*in)>>basename;

		std::string msg_end(" ("+my::tostring(samples_list_.size())+" samples loaded in "+my::tostring(chrono.elapsed())+"s)");
		std::string header(in->get_header());
		header.erase(0,header.find(".. end_of_saved_variables\n",0)+28);
		pso_info_.text(header);
		pso_info_.title("Minimization",'>');
		pso_info_.item(msg+msg_end);
		std::cout<<msg_end<<std::endl;

		delete in;
		in = NULL;
	} else {
		std::string msg("no samples loaded");
		std::cout<<"#"+msg<<std::endl;
		pso_info_.title("Minimization",'>');
		pso_info_.item(msg);

		set_phase_space(P);

		basename = "-" + P.get<std::string>("wf");
		basename+= "-N"  + my::tostring(N);
		basename+= "-m"  + my::tostring(m);
		basename+= "-n"  + my::tostring(n);
		basename+= "-bc" + my::tostring(bc);
		basename+= "-M";
		for(unsigned int i(0);i<M.size();i++){
			basename+= "-"+my::tostring(M(i));
		}
		basename+= "-J";
		for(unsigned int i(0);i<J.size();i++){
			basename+= (J(i)>0?"+":"")+my::tostring(J(i));
		}
	}
	return basename;
}

void VMCMinimization::Minimization::set_phase_space(Parseur& P){
	unsigned int size;
	if(P.find("PS",size)){
		IOFiles load(P.get<std::string>(size),false);
		std::string PS;
		load>>PS;

		std::vector<std::string> ps(my::string_split(PS,'\n'));
		if(ps.size() == Nfreedom_){
			unsigned int k;
			for(unsigned int i(0);i<ps.size();i++){
				ps[i].erase(remove_if(ps[i].begin(),ps[i].end(),isspace),ps[i].end());
				std::vector<std::string> u(my::string_split(ps[i],'U'));
				Vector<double>* vec(new Vector<double>[u.size()]);
				size = 0;
				k = 0;
				for(unsigned int j(0);j<u.size();j++){
					u[j].erase(remove(u[j].begin(),u[j].end(),'['),u[j].end());
					u[j].erase(remove(u[j].begin(),u[j].end(),']'),u[j].end());
					std::vector<std::string> v(my::string_split(u[j],':'));
					if(v.size()==3){
						vec[k]= Vector<double>(my::string2type<double>(v[0]), my::string2type<double>(v[2]), my::string2type<double>(v[1]));
						size += vec[k].size();
						k++;
					} else {
						std::cerr<<"void VMCMinimization::Minimization::set_phase_space(Parseur& P) : each range must be given this way [min:dx:max]"<<std::endl;
					}
				}
				ps_size_ *= size;
				ps_[i].set(size);

				k = 0;
				for(unsigned int j(0);j<u.size();j++){
					for(unsigned int l(0);l<vec[j].size();l++){
						ps_[i](k++) = vec[j](l);
					}
				}
				delete[] vec;
			}

			std::string msg("phase space with "+my::tostring(ps_size_)+" elements");
			std::cout<<"#"+msg<<std::endl;
			pso_info_.item(msg);
			pso_info_.nl();
			pso_info_.lineblock(PS);
		} else {
			std::cerr<<"void VMCMinimization::Minimization::set_phase_space(Parseur& P) : provide Nfreedom_ ranges and remove any blank space and EOL at the EOF"<<std::endl;
		}
	} else {
		std::cerr<<"void VMCMinimization::Minimization::set_phase_space(Parseur& P) : need to provide a file containing the phase space"<<std::endl;
	}
}

VMCMinimization::Minimization::~Minimization(){
	std::cerr<<pso_info_.get()<<std::endl;
	if(ps_){ delete[] ps_; }
	if(s_){ delete s_; }
}

//bool VMCMinimization::Minimization::within_limit(Vector<double> const& x){
//for(unsigned int i(0);i<Nfreedom_;i++){
//if(x(i)<ps_[i](0) || x(i)>ps_[i].back()) { return false; }
//}
//return true;
//}

bool VMCMinimization::Minimization::within_limit(Vector<double> const& x){
	unsigned int i(0);
	unsigned int j(0);
	do{
		if(j==ps_[i].size()){ return false; }
		if(my::are_equal(x(i),ps_[i](j))){ i++; j=0; }
		else { j++; }
	} while (i<Nfreedom_);
	return true;
}

void VMCMinimization::Minimization::save(IOFiles& out) const {
	s_->save(out);

	out.write("Nfreedom",Nfreedom_);
	for(unsigned int i(0);i<Nfreedom_;i++){ out<<ps_[i]; }
	out.write("# samples",samples_list_.size());
	out.add_header()->nl();
	out.add_header()->comment("end_of_saved_variables");
	out.add_header()->text(pso_info_.get());

	samples_list_.set_target();
	while(samples_list_.target_next()){ samples_list_.get().write(out); }
}
/*}*/
