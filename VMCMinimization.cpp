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
	std::cout<<RST::hashtag_line_<<std::endl;
	std::cout<<"#creating VMCMinimization"<<std::endl;
	m_->set(P,path_,basename_);
	if(m_->s_){
		if(m_->s_->get_status() != 4 || P.locked()){
			std::cerr<<__PRETTY_FUNCTION__<<" : something went wrong, status="<<m_->s_->get_status()<<std::endl;
			m_.reset();
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : no System created"<<std::endl;
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
{ set_time(); }

VMCMinimization::VMCMinimization(IOFiles& in):
	time_(""),
	path_(""),
	basename_(""),
	prefix_(""),
	out_(NULL),
	m_(std::make_shared<Minimization>())
{ m_->load(in,path_,basename_); }
/*}*/

/*{public methods*/
void VMCMinimization::refine(){
	std::cout<<RST::hashtag_line_<<std::endl;
	m_->tmax_ = 1;
	double E;
	double dE(0.005);
	while(dE>5e-5){
		E=0;
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			if(m_->samples_.get().get_energy().get_x()<E){
				E = m_->samples_.get().get_energy().get_x();
			}
		}
		E += dE;
		refine(E,dE);
		dE /= 2;
		if(dE<1e-3){ m_->tmax_ *= 2; }
	}
}

void VMCMinimization::refine(double const& E, double const& dE){
	if(m_->samples_.size() && m_->tmax_){
		List<MCSim> best;
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			if(m_->samples_.get().get_energy().get_x()<E){
				best.add_end(m_->samples_.get_ptr());
			}
		}
		total_eval_ = best.size();

		unsigned int maxiter(10);
		std::string msg("refines "+my::tostring(total_eval_)+" samples (max time "+my::tostring(total_eval_*m_->tmax_*maxiter)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E="+my::tostring(dE));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(total_eval_>5 && total_eval_<1000){
			best.set_target();
			progress_ = 0;
			while(best.target_next()){ evaluate_until_precision(best.get().get_param(),false,dE,maxiter); }
			save();
		} else {
			if(total_eval_<1000){ msg = "not enough samples to be usefull, skip the evaluation"; }
			else { msg = "too many samples, would take too much time, skip the evaluation"; }
			std::cout<<"#"<<msg<<std::endl;
			m_->info_.item(msg);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples or tmax_ = 0"<<std::endl; }
}

void VMCMinimization::refine(unsigned int const& nmin, bool const& set_obs, double const& dE, unsigned int const& maxiter){
	if(m_->samples_.size() && m_->tmax_){
		total_eval_ = std::min(nmin,m_->samples_.size());
		List<MCSim> potential_minima;
		List<MCSim> sorted_samples;
		find_minima(0,0.999,sorted_samples,potential_minima);

		std::string msg("refines "+my::tostring(total_eval_)+" samples (max time "+my::tostring(total_eval_*m_->tmax_*maxiter)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E="+my::tostring(dE));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		sorted_samples.set_target();
		progress_ = 0;
		while(sorted_samples.target_next() && progress_<total_eval_){ evaluate_until_precision(sorted_samples.get().get_param(),set_obs,dE,maxiter); }
		save();
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples or tmax_ = 0"<<std::endl; }
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	std::cout<<RST::hashtag_line_<<std::endl;
	std::string msg("complete_analysis called with convergence_criterion="+my::tostring(convergence_criterion));
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);

	m_->samples_.set_target();
	while(m_->samples_.target_next()){
		m_->samples_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	set_time();
	IOFiles out(path_+get_filename()+".jdbin",true);
	m_->save(out,true);
	out<<path_<<basename_;
}

double VMCMinimization::find_minima(unsigned int const& max_pm, double const& range, List<MCSim>& sorted_samples, List<MCSim>& potential_minima) const {
	double E_range(666);
	if(m_->samples_.size()){
		sorted_samples.set();

		/*!finds the MCSim with the minimal energy*/
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			E_range = m_->samples_.get().get_energy().get_x()<E_range?m_->samples_.get().get_energy().get_x():E_range;
		}
		E_range *= range;

		/*!sort by energy all MCSim with energy lower than a threshold*/
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			if(m_->samples_.get().get_energy().get_x()<E_range){
				sorted_samples.add_sort(m_->samples_.get_ptr(),MCSim::sort_by_E);
			}
		}

		if(max_pm){
			/*!select only some local minima*/
			Vector<double>* param(new Vector<double>[max_pm]);

			unsigned int iter;
			bool keep;
			double d_lim(0.8);
			List<MCSim> tmp;
			do{
				d_lim *= d_lim;
				iter = 0;
				tmp.set();
				sorted_samples.set_target();
				while(sorted_samples.target_next() && iter<max_pm){
					param[iter] = sorted_samples.get().get_param();
					keep = true;
					for(unsigned int i(0);i<iter;i++){
						if((param[i]-param[iter]).variance()<d_lim){
							i=iter;
							keep = false;
						}
					}

					if(keep){
						if(!tmp.find_in_sorted_list(sorted_samples.get_ptr(),MCSim::sort_for_merge)){
							tmp.add_after_target(sorted_samples.get_ptr());
							iter++;
						}
					}
				}
			} while( 1.2*iter<max_pm && d_lim>0.01 );

			iter = 0;
			sorted_samples.set_target();
			while(sorted_samples.target_next() && iter++<max_pm){
				if(!tmp.find_in_sorted_list(sorted_samples.get_ptr(),MCSim::sort_for_merge)){
					tmp.add_after_target(sorted_samples.get_ptr());
				}
			}

			potential_minima.set();
			tmp.set_target();
			while(tmp.target_next()){ potential_minima.add_sort(tmp.get_ptr(),MCSim::sort_by_E); }

			delete[] param;
		}
	}
	return E_range;
}

void VMCMinimization::find_and_run_minima(unsigned int const& max_pm, bool const& set_obs, double const& dE){
	if(m_->samples_.size() && m_->tmax_){
		List<MCSim> potential_minima;
		List<MCSim> sorted_samples;
		find_minima(max_pm,0.999,sorted_samples,potential_minima);

		unsigned int maxiter(1);
		total_eval_ = potential_minima.size();
		std::cout<<RST::hashtag_line_<<std::endl;
		std::string msg("compute all ("+my::tostring(m_->obs_.size())+") observables for "+my::tostring(total_eval_)+" minima (max time "+my::tostring(total_eval_*m_->tmax_*maxiter)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		potential_minima.set_target();
		progress_ = 0;
		while(potential_minima.target_next()){ evaluate_until_precision(potential_minima.get().get_param(),set_obs,dE,maxiter); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples or tmax_ = 0"<<std::endl; }
}

void VMCMinimization::find_save_and_plot_minima(unsigned int const& max_pm, IOFiles& w, std::string path, std::string filename) const {
	/*acually it computes r^2 and not r...*/
	if(m_->samples_.size()){
		if(path==""){ path = path_; }
		if(filename==""){ filename =  get_filename(); }

		List<MCSim> potential_minima;
		List<MCSim> sorted_samples;
		double E_range(find_minima(max_pm,0.999,sorted_samples,potential_minima));

		w.write("number of samples",potential_minima.size());

		potential_minima.set_target();
		potential_minima.target_next();
		Vector<double> best_param(potential_minima.get().get_param());
		IOFiles data_Er(path+filename+"-Er.dat",true);
		do{
			data_Er<<(best_param-potential_minima.get().get_param()).norm_squared()<<" "<<potential_minima.get().get_energy().get_x()<<IOFiles::endl;
			potential_minima.get().save(w);
		} while(potential_minima.target_next());

		IOFiles data(path+filename+".dat",true);
		m_->samples_.target_next();
		while(m_->samples_.target_next()){
			data<<m_->samples_.get().get_param()<<" "<<(best_param-m_->samples_.get().get_param()).norm_squared()<<" "<<m_->samples_.get().get_energy()<<IOFiles::endl;
		}

		Gnuplot gp(path,filename);
		gp+="E_range="+my::tostring(E_range);
		gp.multiplot();
		gp.range("x","[:E_range] writeback");
		gp.margin("0.1","0.9","0.5","0.1");
		gp+="plot '"+filename+".dat'    u "+my::tostring(m_->dof_+2)+":"+my::tostring(m_->dof_+1)+":"+my::tostring(m_->dof_+3)+" w xe notitle,\\";
		gp+="     '"+filename+"-Er.dat' u 2:1 lc 4 ps 2 pt 7 t 'selected minima'";
		gp.margin("0.1","0.9","0.9","0.5");
		gp.range("x","restore");
		gp.tics("x");
		gp.key("left Left");
		for(unsigned int i(0);i<m_->dof_;i++){
			gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(m_->dof_+2)+":"+my::tostring(i+1)+":"+my::tostring(m_->dof_+3)+" w xe t '$"+my::tostring(i)+"$'"+(i==m_->dof_-1?"":",\\");
		}
		gp.save_file();
		gp.create_image(true,true);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples"<<std::endl; }
}

void VMCMinimization::explore_around_minima(unsigned int const& max_pm, bool const& set_obs, double const& dE, double const& dx){
	if(m_->samples_.size() && m_->tmax_){
		/*!find the minima and sort by energy*/
		List<MCSim> sorted_samples;
		List<MCSim> potential_minima;
		find_minima(max_pm,0.999,sorted_samples,potential_minima);

		Vector<double> p;
		List<Vector<double> > param;
		std::shared_ptr<Vector<double> > p_ptr;

		/*!select parameters that are close to potential_minima*/
		potential_minima.set_target();
		while(potential_minima.target_next()){
			p = potential_minima.get().get_param();
			for(unsigned int j(0);j<m_->dof_;j++){
				p(j) += dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_in_sorted_list(p_ptr,MCSim::sort_by_param_for_merge)){ param.add_after_target(p_ptr); }

				p(j) -= 2*dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_in_sorted_list(p_ptr,MCSim::sort_by_param_for_merge)){ param.add_after_target(p_ptr); }

				p(j) += dx;
				p_ptr = std::make_shared<Vector<double> >(p);
				if(!param.find_in_sorted_list(p_ptr,MCSim::sort_by_param_for_merge)){ param.add_after_target(p_ptr); }
			}
		}

		unsigned int maxiter(10);
		total_eval_ = param.size();
		std::cout<<RST::hashtag_line_<<std::endl;
		std::string msg("measures "+my::tostring(total_eval_)+" samples close to potential minimas (max time "+my::tostring(m_->tmax_*maxiter*total_eval_)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = "compute ("+my::tostring(m_->obs_.size())+") observables for each samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		param.set_target();
		progress_ = 0;
		while(param.target_next()){ evaluate_until_precision(param.get(),set_obs,dE,maxiter); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples or tmax_ = 0"<<std::endl; }
}

void VMCMinimization::check(unsigned int const& max_pm){
	if(m_->samples_.size()){
		/*!find the minima and sort by energy*/
		List<MCSim> sorted_samples;
		List<MCSim> potential_minima;
		find_minima(max_pm,0.9,sorted_samples,potential_minima);

		Vector<double> param;
		unsigned int i(0);
		sorted_samples.set_target();
		std::cout<<m_->samples_.size()<<std::endl;
		while(sorted_samples.target_next() && i++<max_pm){
			param = sorted_samples.get().get_param();
			std::cout
				<<my::sign(param(0)*param(1)*param(3)*param(4))+1<<" "
				<<my::sign(param(2)*param(3)*param(1)*param(6))+1<<" "
				<<my::sign(param(0)*param(4)*param(5)*param(7))+1<<" "
				<<my::sign(param(2)*param(7)*param(5)*param(6))+1<<std::endl;
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples"<<std::endl; }
}

void VMCMinimization::improve_bad_samples(double const& dE){
	if(m_->samples_.size() && m_->tmax_){
		complete_analysis(1e-5);
		double E(666);
		double tmp;
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			tmp = m_->samples_.get().get_energy().get_x();
			if(tmp<E){ E=tmp; }
		}

		std::vector<Vector<double> > to_improve;
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			if(m_->samples_.get().get_energy().get_x() - 10*m_->samples_.get().get_energy().get_dx() < E && m_->samples_.get().get_energy().get_dx()>dE){
				to_improve.push_back(m_->samples_.get().get_param());
			}
		}

		std::cout<<RST::hashtag_line_<<std::endl;
		std::string msg("improve "+my::tostring(to_improve.size())+" samples (max time "+my::tostring(m_->tmax_*to_improve.size())+"s) with "+RST::math("\\mathrm{d}E>"+my::tostring(dE)));
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		unsigned int nmin(to_improve.size());
#pragma omp parallel for
		for(unsigned int i=0;i<nmin;i++){
			if(omp_get_thread_num()==0){ std::cout<<i<<"/"<<nmin<<std::endl; }
			evaluate(to_improve[i],0);
		}
		complete_analysis(1e-5);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples or tmax_ = 0"<<std::endl; }
}

void VMCMinimization::save_parameters(unsigned int nbest) const {
	List<MCSim> potential_minima;
	List<MCSim> sorted_samples;
	find_minima(0,0.999,sorted_samples,potential_minima);
	nbest = (sorted_samples.size()>nbest+1?nbest:sorted_samples.size()-1);

	set_time();
	IOFiles out(path_+get_filename()+"-parameters.jdbin",true);
	m_->save(out,false);
	out<<path_<<basename_;

	out.write("number of parameters",nbest);

	unsigned int i(0);
	sorted_samples.set_target();
	while(sorted_samples.target_next() && i++<nbest){ out<<sorted_samples.get().get_param(); }
	std::string note(RST::textbf("Maximal Energy :"));
	note += RST::math("E="+my::tostring(sorted_samples.get().get_energy().get_x())) + " ";
	note += RST::math("\\pm\\mathrm{d}E="+my::tostring(sorted_samples.get().get_energy().get_dx()));
	out.add_header()->np();
	out.add_header()->text(note);
}

void VMCMinimization::run_parameters(Parseur& P){
	if(m_->tmax_){
		IOFiles in(P.get<std::string>("param"),false);
		Minimization tmp;
		std::string tmp_path;
		std::string tmp_basename;
		tmp.load(in,tmp_path,tmp_basename);
		Vector<double> param;

		total_eval_ = in.read<unsigned int>();
		progress_ = 0;
		for(unsigned int i(0);i<total_eval_;i++){
			in>>param;
			evaluate_until_precision(param,P.get<bool>("set_obs"),P.get<double>("dE"),P.get<unsigned int>("maxiter"));
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

void VMCMinimization::clean(){
	m_->samples_.set_target();
	while(m_->samples_.target_next()){ m_->samples_.get().set_obs(1); }
}
/*}*/

/*{protected methods*/
std::shared_ptr<MCSim> VMCMinimization::evaluate(Vector<double> const& param, bool const& set_obs){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	List<MCSim>::Node* sample(NULL);
	if(m_->samples_.find_in_sorted_list(sim,sample,MCSim::sort_for_merge)){ sim->copy_S(sample->get()); }
	else {/*!create a new sample*/
		sim->create_S(m_->s_);
		sample = NULL;/*!reset because its value may have been changed*/
	}

	if(sim->is_created()){
		if(set_obs){ for(unsigned int i(1);i<m_->obs_.size();i++){ sim->set_obs(m_->obs_[i]); } }
		sim->run(sample?10:1e6,m_->tmax_);

		if(!sample){
			sim->free_memory();
			/*!if the sample wasn't found before, search again because another
			 * thread may have created it*/
#pragma omp critical(List__global)
			if(!m_->samples_.find_in_sorted_list(sim,sample,MCSim::sort_for_merge)){
				/*!if it is still not found, add the sample to the list*/
				m_->samples_.set_target(sample);
				m_->samples_.add_after_target(sim);
				m_->samples_.set_target();
				sample = NULL;
			}
		}
		if(sample){/*!if the sample exists, merge it with this new measure*/
#pragma omp critical(System__merge)
			sample->get()->merge(sim);
			sim = sample->get();
		}
		return sim;
	} else { return NULL; }
}

void VMCMinimization::evaluate_until_precision(Vector<double> const& param, bool const& set_obs, double const& dE, unsigned int const& maxiter){
	std::shared_ptr<MCSim> sim(NULL);
	unsigned int iter(0);
	std::cout<<++progress_<<"/"<<total_eval_<<" : param : "<<param<<std::endl;
	do {
#pragma omp parallel
		{ sim = evaluate(param,set_obs); }
	} while ( sim.get() && ++iter<maxiter && ( !sim->check_conv(1e-5) || sim->get_energy().get_dx()>dE ) );
	if(sim.get()){
		sim->complete_analysis(1e-5);
		sim->print(0);
	} else { std::cerr<<__PRETTY_FUNCTION__<<std::cout<<" : failed"<<std::endl; }
}
/*}*/
/*}*/

/*{Minimization*/
VMCMinimization::Minimization::~Minimization(){
	if(ps_){ delete[] ps_; }
	if(s_){ delete s_; }
}

void VMCMinimization::Minimization::set(Parseur& P, std::string& path, std::string& basename){
	unsigned int i(0);
	if(P.find("tmax",i,false)){ tmax_ = P.get<unsigned int>(i); }
	if(P.find("load",i,false)){
		Time chrono;
		std::string filename(P.get<std::string>(i));
		std::string msg("loading samples from "+filename);
		std::cout<<"#"+msg<<std::endl;

		IOFiles in(filename,false);
		load(in,path,basename);

		info_.title("Minimization",'>');
		std::string msg_end("("+my::tostring(samples_.size())+" samples loaded in "+my::tostring(chrono.elapsed())+"s)");
		info_.item(msg+msg_end);
		std::cout<<"#"+msg_end<<std::endl;
	} else { create(P,path,basename); }
	if(s_){ s_->print(1); }
}

void VMCMinimization::Minimization::create(Parseur& P, std::string& path, std::string& basename){
	dof_= P.get<unsigned int>("dof");
	if(set_phase_space(P)){
		unsigned int i;
		if(!P.find("M",i,false)){
			std::vector<unsigned int> M(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"));
			P.set("M",M);
		}

		s_  = new System(P);
		Vector<double> tmp(dof_,1.0);

		CreateSystem cs(s_);
		cs.init(&tmp,NULL);
		if(cs.get_status()>3){
			std::cerr<<__PRETTY_FUNCTION__<<" delete s_, status_="<<cs.get_status()<<std::endl;
			delete s_;
			s_ = NULL;
		} else {
			/*!Create all observable and then, set obs_ so it contains every
			 * observable available and set s_ so it contains the minimal
			 * information over the links to compute the energy*/
			cs.create_obs(0);
			obs_ = cs.get_obs();
			s_->set_obs(obs_[0]);
			/*!Sets the bond energies*/
			s_->set_J(cs.get_GenericSystem());

			std::string msg("no samples loaded");
			std::cout<<"#"+msg<<std::endl;
			info_.title("Minimization",'>');
			info_.item(msg);

			path = cs.get_path();
			path+= my::tostring(dof_)+"dof/";
			basename = "-" + cs.get_filename();
			Linux().mkpath(path.c_str());
		}
	}
}

void VMCMinimization::Minimization::load(IOFiles& in, std::string& path, std::string& basename){
	s_ = new System(in);
	in>>dof_;
	Vector<double> tmp(dof_,1.0);

	/*!Get rid of all observables so if new ones have been implemented, they
	 * will be available via the use of CreateSystem. Once they are all
	 * created, set obs_ so it contains every observable available and set s_
	 * so it contains the minimal information over the links to compute the
	 * energy*/
	s_->set_obs(0);
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.create_obs(0);
	obs_ = cs.get_obs();
	s_->set_obs(obs_[0]);

	ps_size_ = 1;
	ps_= new Vector<double>[dof_];
	for(unsigned int i(0);i<dof_;i++){
		ps_[i] = in.read<Vector<double> >();
		ps_size_ *= ps_[i].size();
	}

	unsigned int n_samples(in.read<unsigned int>());
	unsigned int iter(0);
	std::string msg("loading progress : ");
	while(iter++<n_samples){
		samples_.add_end(std::make_shared<MCSim>(in));
		if(!(iter%1000)){ my::display_progress(iter,n_samples,msg); }
	}
	in>>path>>basename;

	std::string header(in.get_header());
	header.erase(0,header.find(".. end_of_saved_variables\n",0)+27);
	info_.text(header);
}

bool VMCMinimization::Minimization::set_phase_space(Parseur const& P){
	unsigned int size;
	if(P.find("PS",size,true)){
		IOFiles load(P.get<std::string>(size),false);
		std::string PS;
		load>>PS;

		ps_size_ = 1;
		std::vector<std::string> ps(my::string_split(PS,'\n'));
		if(ps.size() == dof_){
			unsigned int k;

			if(ps_){ delete[] ps_; }
			ps_ = new Vector<double>[dof_];
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
			return true;
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : provide "<<dof_<<" ranges and remove any blank space and EOL at the EOF"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : need to provide a file containing the phase space"<<std::endl; }
	return false;
}

//bool VMCMinimization::Minimization::within_limit(Vector<double> const& x){
//for(unsigned int i(0);i<dof_;i++){
//if(x(i)<ps_[i](0) || x(i)>ps_[i].back()) { return false; }
//}
//return true;
//}

bool VMCMinimization::Minimization::within_limit(Vector<double> const& x) const {
	unsigned int i(0);
	unsigned int j(0);
	do{
		if(j==ps_[i].size()){ return false; }
		if(my::are_equal(x(i),ps_[i](j))){ i++; j=0; }
		else { j++; }
	} while (i<dof_);
	return true;
}

void VMCMinimization::Minimization::save(IOFiles& out, bool const& all) const {
	/*!saves a system that has not been measured but it is required for the
	 * eventual call of System(IOFiles& r).*/
	s_->save(out);

	out.write("dof",dof_);
	for(unsigned int i(0);i<dof_;i++){ out<<ps_[i]; }
	out.write("# samples",samples_.size());

	Vector<double> tmp(dof_,1.0);
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	out.add_header()->np();
	out.add_header()->text(cs.get_system_info().get());
	out.add_header()->np();

	if(all){
		double E(0);
		List<MCSim>::Node* best(NULL);
		samples_.set_target();
		while(samples_.target_next()){
			samples_.get().write(out);
			if(samples_.get().get_energy().get_x()<E){
				E = samples_.get().get_energy().get_x();
				best = samples_.get_target();
			}
		}
		std::string p("(");
		for(unsigned int i(0);i<best->get()->get_param().size()-1;i++){ p += my::tostring(best->get()->get_param()(i))+","; }
		p += my::tostring(best->get()->get_param().back())+")";
		out.add_header()->def("Best parameter","p="+p);
		out.add_header()->def("Best energy","E="+my::tostring(best->get()->get_energy().get_x())+", dE="+my::tostring(best->get()->get_energy().get_dx()));
		out.add_header()->np();
	}

	out.add_header()->comment("end_of_saved_variables");
	out.add_header()->text(info_.get());
}
/*}*/
