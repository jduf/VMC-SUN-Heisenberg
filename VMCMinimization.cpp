#include "VMCMinimization.hpp"

/*{VMCMinimization*/
/*{construcor*/
VMCMinimization::VMCMinimization(Parseur& P):
	time_(""),
	path_(""),
	basename_(""),
	prefix_("MIN"),
	out_(NULL),
	m_(std::make_shared<Minimization>())
{
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#creating VMCMinimization"<<std::endl;
	m_->set(P,path_,basename_);
	if(m_->s_->get_status() != 3 || P.locked()){
		std::cerr<<__PRETTY_FUNCTION__<<" : something went wrong, status="<<m_->s_->get_status()<<std::endl;
		m_.reset();
	}
}

VMCMinimization::VMCMinimization(VMCMinimization const& vmcm, std::string const& prefix):
	time_(""),
	path_(vmcm.path_),
	basename_(vmcm.basename_),
	prefix_(prefix),
	out_(NULL),
	m_(vmcm.m_)
{}

VMCMinimization::VMCMinimization(IOFiles& in):
	time_(""),
	path_(""),
	basename_(""),
	prefix_(""),
	out_(NULL),
	m_(std::make_shared<Minimization>())
{
	m_->load(in,path_,basename_);
}
/*}*/

/*{public methods*/
void VMCMinimization::refine(){
	std::cout<<"#######################"<<std::endl;
	m_->tmax_ = 5;
	double E;
	double dE(0.1);
	while(dE>5e-5){
		E=0;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_MCS()->get_energy().get_x()<E){
				E = m_->samples_list_.get().get_MCS()->get_energy().get_x();
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
		unsigned int maxiter(10);
		List<MCSim> best;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_MCS()->get_energy().get_x()<E){
				best.add_end(m_->samples_list_.get_ptr());
			}
		}
		unsigned int N(best.size());

		std::string msg("refines "+my::tostring(N)+" samples (max time "+my::tostring(N*m_->tmax_ * maxiter)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E="+my::tostring(dE));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(N>5 && N<700){
			best.set_target();
			while(best.target_next()){ evaluate_until_precision(best.get().get_param(),dE,-1,maxiter); }
		} else {
			if(N<700){ msg = "not enough data to be usefull, skip the evaluation"; }
			else { msg = "too many data, would take too much time, skip the evaluation"; }
			std::cout<<"#"<<msg<<std::endl;
			m_->info_.item(msg);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl; }
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	std::cout<<"#######################"<<std::endl;
	std::string msg("complete_analysis called with convergence_criterion="+my::tostring(convergence_criterion));
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);

	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		m_->samples_list_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	set_time();
	IOFiles out(path_+get_filename()+".jdbin",true);
	m_->save(out);
	out<<path_<<basename_;
}

void VMCMinimization::find_minima(unsigned int const& max_n_minima, List<MCSim>& list_min, Vector<double>& best_param, double& E_range) const {
	Vector<double>* param(new Vector<double>[max_n_minima]);

	/*!finds the MCSim with the minimal energy*/
	double tmp;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		tmp = m_->samples_list_.get().get_MCS()->get_energy().get_x();
		if(tmp<E_range){ E_range=tmp; }
	}
	E_range *= 0.99;

	List<MCSim> sorted_list;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		if(m_->samples_list_.get().get_MCS()->get_energy().get_x()<E_range){
			sorted_list.add_sort(m_->samples_list_.get_ptr(),MCSim::sort_by_E); 
		}
	}

	unsigned int local_min;
	bool keep;
	double d_lim(0.8);
	do{
		d_lim *= d_lim;
		local_min = 0;
		list_min.set();
		sorted_list.set_target();
		while(sorted_list.target_next() && local_min<max_n_minima){
			param[local_min] = sorted_list.get().get_param();
			keep = true;
			for(unsigned int i(0);i<local_min;i++){
				if((param[i]-param[local_min]).variance()<d_lim){
					i=local_min;
					keep = false;
				}
			}

			if(keep){
				list_min.add_end(sorted_list.get_ptr());
				local_min++;
			}
		}
	} while( 1.2*local_min<max_n_minima && d_lim>0.01 );

	best_param = param[0];
	delete[] param;
}

void VMCMinimization::find_save_and_plot_minima(unsigned int const& max_n_minima, IOFiles& w, std::string path, std::string filename) const {
	/*acually it computes r^2 and not r...*/
	if(m_->samples_list_.size()){
		if(path==""){ path = path_; }
		if(filename==""){ filename =  get_filename(); }

		List<MCSim> list_min;
		Vector<double> best_param;
		double E_range;
		find_minima(max_n_minima,list_min,best_param,E_range);

		IOFiles data(path+filename+".dat",true);
		IOFiles data_Er(path+filename+"-Er.dat",true);
		w.write("number of minima",list_min.size());

		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			data<<m_->samples_list_.get().get_param()<<" "<<(best_param-m_->samples_list_.get().get_param()).norm_squared()<<" "<<m_->samples_list_.get().get_MCS()->get_energy()<<IOFiles::endl;
		}

		list_min.set_target();
		while(list_min.target_next()){
			data_Er<<(best_param-list_min.get().get_param()).norm_squared()<<" "<<list_min.get().get_MCS()->get_energy().get_x()<<IOFiles::endl;
			list_min.get().save(w);
		}

		Gnuplot gp(path,filename);
		gp+="E_range="+my::tostring(E_range);
		gp.multiplot();
		gp.range("x","[:E_range] writeback");
		gp.margin("0.1","0.9","0.5","0.10");
		gp+="plot '"+filename+".dat'        u "+my::tostring(m_->dof_+2)+":"+my::tostring(m_->dof_+1)+":"+my::tostring(m_->dof_+3)+" w xe           notitle,\\";
		gp+="     '"+filename+"-Er.dat'     u 2:1   lc 4 ps 2 pt 7 t 'selected minima'";
		gp.margin("0.1","0.9","0.9","0.50");
		gp.tics("x");
		gp.range("x","restore");
		gp.key("left Left");
		for(unsigned int i(0);i<m_->dof_;i++){
			gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(m_->dof_+2)+":"+my::tostring(i+1)+":"+my::tostring(m_->dof_+3)+" w xe t '$"+my::tostring(i)+"$'"+(i==m_->dof_-1?"":",\\");
		}
		gp.save_file();
		gp.create_image(true,true);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl; }
}

void VMCMinimization::explore_around_minima(unsigned int const& max_n_minima, int const& nobs, double const& dE, double const& dx){
	if(m_->samples_list_.size()){
		unsigned int maxiter(10);
		List<MCSim> sorted_list;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			sorted_list.add_sort(m_->samples_list_.get_ptr(),MCSim::sort_by_E);
		}

		auto sort_by_param = [](Vector<double> const& a, Vector<double> const& b){
			for(unsigned int i(0);i<a.size();i++){
				if(a(i) - b(i) > 0.0001){ return 0; }
				if(a(i) - b(i) <-0.0001){ return 1; }
			}
			return 2;
		};

		List<Vector<double> > param;
		Vector<double> p;
		std::shared_ptr<Vector<double> > p_ptr;
		unsigned int i(0);
		sorted_list.set_target();
		while(sorted_list.target_next() && i++<max_n_minima){
			p = sorted_list.get().get_param();
			for(unsigned int j(0);j<m_->dof_;j++){
				p(j) += dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_sorted(p_ptr,sort_by_param)){ param.add_after_target(p_ptr); }

				p(j) -= 2*dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_sorted(p_ptr,sort_by_param)){ param.add_after_target(p_ptr); }

				p(j) += dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_sorted(p_ptr,sort_by_param)){ param.add_after_target(p_ptr); }
			}
		}

		Vector<double> best_param;
		double E_range;
		sorted_list.set();
		find_minima(max_n_minima,sorted_list,best_param,E_range);

		sorted_list.set_target();
		i=0;
		while(sorted_list.target_next() && i++<max_n_minima){
			p = sorted_list.get().get_param();
			for(unsigned int j(0);j<m_->dof_;j++){
				p(j) += dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_sorted(p_ptr,sort_by_param)){ param.add_after_target(p_ptr); }

				p(j) -= 2*dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_sorted(p_ptr,sort_by_param)){ param.add_after_target(p_ptr); }

				p(j) += dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_sorted(p_ptr,sort_by_param)){ param.add_after_target(p_ptr); }
			}
		}

		std::cout<<"#######################"<<std::endl;
		std::string msg("measures "+my::tostring(param.size())+" samples close to local minimas (max time "+my::tostring(m_->tmax_ * maxiter * param.size())+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = "compute "+my::tostring(nobs)+" observables for each samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E="+my::tostring(dE));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		param.set_target();
		while(param.target_next()){ evaluate_until_precision(param.get(),dE,nobs,maxiter); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl; }
}

void VMCMinimization::improve_bad_samples(double const& dE){
	if(m_->samples_list_.size()){
		unsigned int maxiter(10);
		double E(666);
		double tmp;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			tmp = m_->samples_list_.get().get_MCS()->get_energy().get_x();
			if(tmp<E){ E=tmp; }
		}
		E *= 0.9;

		List<MCSim> to_improve;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_MCS()->get_energy().get_dx()>dE && m_->samples_list_.get().get_MCS()->get_energy().get_x()<E){
				to_improve.add_end(m_->samples_list_.get_ptr());
			}
		}

		std::cout<<"#######################"<<std::endl;
		std::string msg("improve "+my::tostring(to_improve.size())+" samples (max time "+my::tostring(m_->tmax_ * maxiter * to_improve.size())+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E="+my::tostring(dE));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		to_improve.set_target();
		while(to_improve.target_next()){ evaluate_until_precision(to_improve.get().get_param(),dE,0,maxiter); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl; }
}

void VMCMinimization::find_and_run_minima(unsigned int const& max_n_minima, int const& nobs, double const& dE){
	if(m_->samples_list_.size()){
		unsigned int maxiter(10);
		double E_range;
		List<MCSim> list_min;
		Vector<double> best_param;
		find_minima(max_n_minima,list_min,best_param,E_range);

		std::cout<<"#######################"<<std::endl;
		std::string msg("compute "+my::tostring(nobs)+" observables for "+my::tostring(list_min.size())+" local minima (max time "+my::tostring(m_->tmax_ * maxiter)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E="+my::tostring(dE));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		list_min.set_target();
		while(list_min.target_next()){ evaluate_until_precision(list_min.get().get_param(),dE,nobs,maxiter); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl; }
}

void VMCMinimization::print() const {
	std::cout<<"Print whole history ("<< m_->samples_list_.size()<<")"<<std::endl;
	while( m_->samples_list_.target_next() ){
		std::cout<<m_->samples_list_.get_ptr()<<" "<<m_->samples_list_.get().get_param()<<" "<<m_->samples_list_.get().get_MCS()->get_energy()<<std::endl;
	}
}
/*}*/

/*{protected methods*/
std::shared_ptr<MCSim> VMCMinimization::evaluate(Vector<double> const& param, int const& nobs){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	List<MCSim>::Node* exists_sample(NULL);
#pragma omp critical(samples_list_)
	{
		if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
			sim->copy_S(m_->samples_list_.get().get_MCS());
			exists_sample = m_->samples_list_.get_target();
		} else { sim->create_S(m_->s_); }
		m_->samples_list_.set_target();
	}
	if(sim->is_created()){
		sim->set_observables(m_->obs_,nobs);
		sim->run(exists_sample?10:1e6,m_->tmax_);
#pragma omp critical(samples_list_)
		{
			if(exists_sample){
				m_->samples_list_.set_target(exists_sample);
				m_->samples_list_.handle_twin(sim,MCSim::merge);
				//sim = exists_sample->get();/*should be equivalent*/
				sim = m_->samples_list_.get_ptr();
			} else {
				/*search because the sample may have been created by another
				 *thread*/
				if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
					m_->samples_list_.handle_twin(sim,MCSim::merge);
					sim = m_->samples_list_.get_ptr();
				} else { m_->samples_list_.add_after_target(sim); }
				m_->samples_list_.set_target();
			}
		}
		sim->free_memory();
		return sim;
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : not valid parameter : "<<param<<std::endl;
		return NULL;
	}
}

void VMCMinimization::evaluate_until_precision(Vector<double> const& param, double const& dE, int const& nobs, unsigned int const& maxiter){
	std::shared_ptr<MCSim> sim(NULL);
	unsigned int iter(0);
	do {
		if(sim.get()){ std::cout<<iter<<" iter "<<param<<" "<<sim->get_MCS()->get_energy()<<std::endl; }
#pragma omp parallel
		{ sim = evaluate(param,nobs); }
	} while ( sim.get() && ++iter<maxiter && ( !sim->check_conv(1e-5) || sim->get_MCS()->get_energy().get_dx()>dE ) );
	if(sim.get()){ sim->complete_analysis(1e-5); }
}
/*}*/
/*}*/

/*{Minimization*/
void VMCMinimization::Minimization::set(Parseur& P, std::string& path, std::string& basename){
	unsigned int i(0);
	if(P.find("load",i,false)){
		Time chrono;
		std::string filename(P.get<std::string>(i));
		std::string msg("loading samples from "+filename);
		std::cout<<"#"+msg<<std::flush;

		IOFiles in(filename,false);
		std::string n_samples(load(in,path,basename));

		info_.title("Minimization",'>');
		std::string msg_end(" ("+n_samples+" samples loaded in "+my::tostring(chrono.elapsed())+"s)");
		info_.item(msg+msg_end);
		std::cout<<msg_end<<std::endl;
	} else { create(P,path,basename); }

	if(P.find("tmax",i,false)){ tmax_ = P.get<unsigned int>("tmax"); }
	else {
		tmax_ = 1;
		std::string msg("assume "+RST::math("t_{max} = "+my::tostring(tmax_)+"s"));
		std::cout<<"#"+msg<<std::endl;
		info_.item(msg);
	}
}

void VMCMinimization::Minimization::create(Parseur& P, std::string& path, std::string& basename){
	unsigned int i;
	if(!P.find("M",i,false)){
		std::vector<unsigned int> M(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"));
		P.set("M",M);
	}
	s_  = new System(P);
	dof_= P.get<unsigned int>("dof");
	ps_ = new Vector<double>[dof_];

	/*!Sets obs_ and gives s_ the list of nearest neighbour links*/
	Vector<double> tmp(dof_,1.0);
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.set_observables(-1);
	obs_ = cs.get_GS()->get_obs();
	s_->set_observables(obs_,0);

	std::string msg("no samples loaded");
	std::cout<<"#"+msg<<std::endl;
	info_.title("Minimization",'>');
	info_.item(msg);

	set_phase_space(P);

	Linux command;
	path = cs.get_path();
	path+= my::tostring(dof_)+"dof/";
	basename = "-" + cs.get_filename();
	command.mkpath(path.c_str());
}

std::string VMCMinimization::Minimization::load(IOFiles& in, std::string& path, std::string& basename){
	s_ = new System(in);
	in>>dof_;
	ps_= new Vector<double>[dof_];

	/*!Sets obs_ (s_ should already know the list of nearest neighbour links)*/
	Vector<double> tmp(dof_,1.0);
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.set_observables(-1);
	obs_ = cs.get_GS()->get_obs();

	ps_size_ = 1;
	for(unsigned int i(0);i<dof_;i++){
		ps_[i] = in.read<Vector<double> >();
		ps_size_ *= ps_[i].size();
	}

	unsigned int n_samples(in.read<unsigned int>());
	while(n_samples--){ samples_list_.add_end(std::make_shared<MCSim>(in)); }
	in>>path>>basename;

	std::string header(in.get_header());
	header.erase(0,header.find(".. end_of_saved_variables\n",0)+28);
	info_.text(header);
	return my::tostring(samples_list_.size());
}

void VMCMinimization::Minimization::set_phase_space(Parseur const& P){
	unsigned int size;
	if(P.find("PS",size,true)){
		IOFiles load(P.get<std::string>(size),false);
		std::string PS;
		load>>PS;

		ps_size_ = 1;
		std::vector<std::string> ps(my::string_split(PS,'\n'));
		if(ps.size() == dof_){
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
					} else { std::cerr<<__PRETTY_FUNCTION__<<" : each range must be given this way [min:dx:max]"<<std::endl; }
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

			std::string msg("phase space contains "+my::tostring(ps_size_)+" values");
			std::cout<<"#"+msg<<std::endl;
			info_.item(msg);
			info_.nl();
			info_.lineblock(PS);
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : provide "<<dof_<<" ranges and remove any blank space and EOL at the EOF"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : need to provide a file containing the phase space"<<std::endl; }
}

VMCMinimization::Minimization::~Minimization(){
	if(ps_){ delete[] ps_; }
	if(s_){ delete s_; }
}

//bool VMCMinimization::Minimization::within_limit(Vector<double> const& x){
//for(unsigned int i(0);i<dof_;i++){
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
	} while (i<dof_);
	return true;
}

void VMCMinimization::Minimization::save(IOFiles& out) const {
	/*!Saves a system that has not been measured but it is required for the
	 * eventual call of System(IOFiles& r).*/
	s_->save_input(out);
	s_->save_output(out);

	out.write("dof",dof_);
	for(unsigned int i(0);i<dof_;i++){ out<<ps_[i]; }
	out.write("# samples",samples_list_.size());
	out.add_header()->nl();
	out.add_header()->comment("end_of_saved_variables");
	out.add_header()->text(info_.get());

	samples_list_.set_target();
	while(samples_list_.target_next()){ samples_list_.get().write(out); }
}
/*}*/
