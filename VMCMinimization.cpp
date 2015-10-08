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
		List<MCSim> best;
		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_MCS()->get_energy().get_x()<E){
				best.add_end(m_->samples_list_.get_ptr());
			}
		}
		unsigned int N(best.size());

		std::string msg("refine "+my::tostring(N)+" samples ("+my::tostring(E)+","+my::tostring(dE)+","+my::tostring(m_->tmax_)+")");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(N>5 && N<700){
			unsigned int iter;
			best.set_target();
			while(best.target_next()){
				iter = 0;
				do {
#pragma omp parallel
					{ evaluate(best.get().get_param()); }
				} while( iter++<10 && ( !best.get().check_conv(1e-5) || best.get().get_MCS()->get_energy().get_dx()>dE ) );
			}
		} else {
			if(N<700){ msg = "not enough data to be usefull, skip the evaluation"; }
			else { msg = "too many data, would take too much time, skip the evaluation"; }
			std::cout<<"#"<<msg<<std::endl;
			m_->info_.item(msg);
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl;
	}
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
	double E(666);
	double tmp;

	auto sort_by_r = [&](MCSim const& a, MCSim const& b){
		double tmp_a((a.get_param()-best_param).norm_squared());
		double tmp_b((b.get_param()-best_param).norm_squared());
		if(my::are_equal(tmp_a,tmp_b)){ return 2; }
		if(tmp_a>tmp_b){ return 0; }
		if(tmp_a<tmp_b){ return 1; }
		return 2;
	};

	/*!finds the MCSim with the minimal energy*/
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		tmp = m_->samples_list_.get().get_MCS()->get_energy().get_x();
		if(tmp<E){
			E=tmp;
			best_param = m_->samples_list_.get().get_param();
		}
	}
	E_range = E*0.99;

	/*!lists all MCSim with energy below E_range and keep only one per r*/
	List<MCSim> list_sorted_r;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		E=m_->samples_list_.get().get_MCS()->get_energy().get_x();
		if(E<E_range){
			if(list_sorted_r.find_sorted(m_->samples_list_.get_ptr(),sort_by_r)){
				if(list_sorted_r.get().get_MCS()->get_energy().get_x() > E){ list_sorted_r.get_ptr() = m_->samples_list_.get_ptr(); }
			} else { list_sorted_r.add_after_target(m_->samples_list_.get_ptr()); }
		}
	}

	/*!finds the minima*/
	unsigned int ao(1);
	do{
		ao *= 2;
		list_min.set();

		bool keep;
		List<MCSim> list_tmp;
		std::vector<double> E_tmp;
		for(unsigned int i(0);i<list_sorted_r.size();i++){
			tmp = (best_param-list_sorted_r[i].get_param()).norm_squared();
			E = list_sorted_r.get().get_MCS()->get_energy().get_x();
			list_tmp.add_end(list_sorted_r.get_ptr());
			keep = true;
			for(unsigned int j(i>ao?i-ao:1);j<i+ao && j<list_sorted_r.size();j++){
				if(list_sorted_r[j].get_MCS()->get_energy().get_x()<E){
					j = list_sorted_r.size();
					keep = false;
				}
			}
			if(keep){ E_tmp.push_back(E); }
			else { list_tmp.pop_end(); }
		}

		unsigned int i(1);
		list_tmp.set_target();
		list_tmp.target_next();
		list_min.add_start(list_tmp.get_ptr());
		/*!the condition's order is important because the last entry should
		 * be targeted otherwise the last entry might be forgotten.*/
		while(list_tmp.target_next() && i<E_tmp.size()-1){
			if( E_tmp[i-1]>E_tmp[i] && E_tmp[i]<E_tmp[i+1] ){ list_min.add_sort(list_tmp.get_ptr(),MCSim::sort_by_E); }
			i++;
		}
		if(E_tmp[i-1]>E_tmp[i]){ list_min.add_sort(list_tmp.get_ptr(),MCSim::sort_by_E); }
	} while ( list_min.size()>max_n_minima );
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
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl;
	}
}

void VMCMinimization::find_and_run_minima(unsigned int const& max_n_minima){
	if(m_->samples_list_.size()){
		std::cout<<"#######################"<<std::endl;
		std::string msg("compute correlations and long range correlations for minima");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		List<MCSim> list_min;
		Vector<double> best_param;
		double E_range;
		find_minima(max_n_minima,list_min,best_param,E_range);

		msg = "found "+my::tostring(list_min.size())+" local minima";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		list_min.set_target();
		while(list_min.target_next()){
#pragma omp parallel
			{ evaluate(list_min.get().get_param(),2); }
			list_min.get().complete_analysis(1e-5);
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl;
	}
}

void VMCMinimization::print() const {
	std::cout<<"Print whole history ("<< m_->samples_list_.size()<<")"<<std::endl;
	while( m_->samples_list_.target_next() ){
		std::cout<<m_->samples_list_.get_ptr()<<" "<<m_->samples_list_.get().get_param()<<" "<<m_->samples_list_.get().get_MCS()->get_energy()<<std::endl;
	}
}
/*}*/

/*{protected methods*/
std::shared_ptr<MCSim> VMCMinimization::evaluate(Vector<double> const& param, unsigned int const& which){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	bool tmp_test;
#pragma omp critical(samples_list_)
	{
		if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
			tmp_test = true;
			sim->copy_S(m_->samples_list_.get().get_MCS());
			/*will need to provide a way to set correct observables*/
			//sim->set_observables(which);
		} else {
			tmp_test = false;
			sim->create_S(m_->s_,which);
		}
		m_->samples_list_.set_target();
	}
	if(sim->is_created()){
		sim->run(tmp_test?10:1e6,m_->tmax_);
#pragma omp critical(samples_list_)
		{
			if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
				m_->samples_list_.handle_twin(sim,MCSim::merge);
				sim = m_->samples_list_.get_ptr();
			} else {
				m_->samples_list_.add_after_target(sim);
			}
			m_->samples_list_.set_target();
		}
		sim->free_memory();
		return sim;
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : not valid parameter : "<<param<<std::endl;
		return NULL;
	}
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
	} else {
		create(P,path,basename);
	}

	if(P.find("tmax",i,false)){ tmax_ = P.get<unsigned int>("tmax"); }
	else {
		tmax_ = 1;
		std::string msg("assume tmax=1s");
		std::cout<<"#"+msg<<std::endl;
		info_.item(msg);
	}
}

void VMCMinimization::Minimization::create(Parseur& P, std::string& path, std::string& basename){
	s_  = new System(P);
	dof_= P.get<unsigned int>("dof");
	ps_ = new Vector<double>[dof_];

	/*!the next block is required to configure J correctly so that path and
	 * filename are correct*/
	Vector<double> tmp;
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.set(s_);

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

	/*!the next block is required to compute the J setup*/
	Vector<double> tmp;
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.set(s_);

	ps_size_ = 1;
	for(unsigned int i(0);i<dof_;i++){
		ps_[i] = in.read<Vector<double> >();
		ps_size_ *= ps_[i].size();
	}

	unsigned int n_samples(in.read<int>());
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
					} else {
						std::cerr<<__PRETTY_FUNCTION__<<" : each range must be given this way [min:dx:max]"<<std::endl;
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

			std::string msg("phase space contains "+my::tostring(ps_size_)+" values");
			std::cout<<"#"+msg<<std::endl;
			info_.item(msg);
			info_.nl();
			info_.lineblock(PS);
		} else {
			std::cerr<<__PRETTY_FUNCTION__<<" : provide "<<dof_<<" ranges and remove any blank space and EOL at the EOF"<<std::endl;
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : need to provide a file containing the phase space"<<std::endl;
	}
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
	/*{Description*/
	/*!will save "empty" E_,corr_,lr_corr_ but it is required if one want to
	 * create a System by calling System(IOFiles& r) later. Note that as s_ is
	 * an instance of System, save_input will not save any parameter*/
	/*}*/
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
