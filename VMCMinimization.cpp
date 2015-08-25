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
	set_time();
	if(m_->s_->get_status() != 3 || P.locked()){
		std::cout<<m_->s_->get_status()<<std::endl;
		m_.reset();
		std::cerr<<__PRETTY_FUNCTION__<<" : something went wrong"<<std::endl;
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
	double E;
	double dE(0.1);
	while(dE>5e-5){
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

#pragma omp parallel for schedule(dynamic,10)
		for(unsigned int i=0;i<N;i++){
			std::shared_ptr<MCSim> sim;
#pragma omp critical
			{
				best.target_next();
				sim = std::make_shared<MCSim>(best.get().get_param());
				/*need to do call copy_S because best.get().get_S() has an
				 * undefined Ainv_, EVec_[i>0] ... due to the free_memory()
				 * call*/
				sim->copy_S(best.get().get_S());
			}
			while(!sim->check_conv(1e-5) || sim->get_S()->get_energy().get_dx()>dE) { sim->run(0,m_->tmax_); }
			/*Merge this new evaluation*/
#pragma omp critical(samples_list_)
			{
				m_->samples_list_.set_target();
				if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
					m_->samples_list_.merge_with_target(sim,MCSim::merge);
					m_->samples_list_.get().get_S()->get_energy().complete_analysis(1e-5);
				} else {
					std::cerr<<__PRETTY_FUNCTION__<<" : can't find sim in m_->samples_list_"<<std::endl;
				}
			}
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl;
	}
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	set_time();
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#complete_analysis called with convergence_criterion="<<convergence_criterion<<std::endl;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		m_->samples_list_.get().complete_analysis(convergence_criterion);
	}
}

void VMCMinimization::save() const {
	IOFiles out(path_+get_filename()+".jdbin",true);
	m_->save(out);
	out<<path_<<basename_;
}

void VMCMinimization::find_minima(unsigned int const& max_n_minima, List<MCSim>& list_min, Vector<double>& param, double& E_range, Interpolation<double>* interp_Er) const {
	double E(666);
	double tmp;

	auto sort_by_r = [&](MCSim const& a, MCSim const& b){
		double tmp_a((a.get_param()-param).norm_squared());
		double tmp_b((b.get_param()-param).norm_squared());
		if(my::are_equal(tmp_a,tmp_b)){ return 2; }
		if(tmp_a>tmp_b){ return 0; }
		if(tmp_a<tmp_b){ return 1; }
		return 2;
	};

	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		tmp = m_->samples_list_.get().get_S()->get_energy().get_x();
		if(tmp<E){
			E=tmp;
			param = m_->samples_list_.get().get_param();
		}
	}
	E_range = E*0.99;

	List<MCSim> list_sorted_r;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		E=m_->samples_list_.get().get_S()->get_energy().get_x();
		if(E<E_range){
			if(list_sorted_r.find_sorted(m_->samples_list_.get_ptr(),sort_by_r)){
				if(list_sorted_r.get().get_S()->get_energy().get_x() > E){ list_sorted_r.get_ptr() = m_->samples_list_.get_ptr(); }
			} else { list_sorted_r.add_after_target(m_->samples_list_.get_ptr()); }
		}
	}

	unsigned int ao(1);
	do{
		ao *= 2;
		list_min.set();
		if(interp_Er){ interp_Er->set_data(); }

		bool keep;
		List<MCSim> list_tmp;
		std::vector<double> E_tmp;
		for(unsigned int i(0);i<list_sorted_r.size();i++){
			tmp = (param-list_sorted_r[i].get_param()).norm_squared();
			E = list_sorted_r.get().get_S()->get_energy().get_x();
			list_tmp.add_end(list_sorted_r.get_ptr());
			keep = true;
			for(unsigned int j(i>ao?i-ao:1);j<i+ao && j<list_sorted_r.size();j++){
				if(list_sorted_r[j].get_S()->get_energy().get_x()<E){
					j = list_sorted_r.size();
					keep = false;
				}
			}
			if(keep){
				if(interp_Er){ interp_Er->add_data(tmp,E); }
				E_tmp.push_back(E);
			} else {
				list_tmp.pop_end();
			}
		}
		if(interp_Er){
			double dx(0.3);
			interp_Er->compute_weights(dx,interp_Er->get_N());
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
		Vector<double> param;
		double E_range;

		//will certainly get rid of this interpolation...
		Interpolation<double> interp_Er(1);
		interp_Er.select_basis_function(7);

		find_minima(max_n_minima,list_min,param,E_range,&interp_Er);

		IOFiles data(path+filename+".dat",true);
		IOFiles data_Er(path+filename+"-Er.dat",true);
		IOFiles data_interp(path+filename+"-interp.dat",true);
		w.write("number of minima",list_min.size());

		m_->samples_list_.set_target();
		while(m_->samples_list_.target_next()){
			data<<m_->samples_list_.get().get_param()<<" "<<(param-m_->samples_list_.get().get_param()).norm_squared()<<" "<<m_->samples_list_.get().get_S()->get_energy()<<IOFiles::endl;
		}

		list_min.set_target();
		double tmp;
		double r_max(0);
		while(list_min.target_next()){
			tmp = (param-list_min.get().get_param()).norm_squared();
			data_Er<<tmp<<" "<<list_min.get().get_S()->get_energy().get_x()<<IOFiles::endl;
			list_min.get().save(w);
			if(r_max<tmp){ r_max=tmp; }
		}

		Vector<double> range(0,r_max,0.05);
		for(unsigned int i(0);i<range.size();i++){
			data_interp<<range(i)<<" "<<interp_Er(range(i))<<IOFiles::endl;
		}

		Gnuplot gp(path,filename);
		gp+="E_range="+my::tostring(E_range);
		gp.multiplot();
		gp.range("x","[:E_range] writeback");
		gp.margin("0.1","0.9","0.5","0.10");
		gp+="plot '"+filename+".dat'        u "+my::tostring(m_->dof_+2)+":"+my::tostring(m_->dof_+1)+":"+my::tostring(m_->dof_+3)+" w xe           notitle,\\";
		gp+="     '"+filename+"-interp.dat' u 2:1   w l            notitle,\\";
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
		List<MCSim> list_min;
		Vector<double> param;
		double E_range;

		find_minima(max_n_minima,list_min,param,E_range);
		std::cout<<"#######################"<<std::endl;
		std::string msg("found "+my::tostring(list_min.size())+" local minima");
		std::cout<<"#"<<msg<<std::endl;
		m_->pso_info_.item(msg);

		list_min.set_target();
		while(list_min.target_next()){ 
			evaluate(list_min.get().get_param(),2);
			list_min.get().complete_analysis(1e-5);
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : there is no data"<<std::endl;
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
/*}*/

/*{protected methods*/
std::shared_ptr<MCSim> VMCMinimization::evaluate(Vector<double> const& param, unsigned int const& which){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
	bool tmp_test;
#pragma omp critical(samples_list_)
	{
		if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
			tmp_test = true;
			sim->copy_S(m_->samples_list_.get().get_S());
		} else {
			tmp_test = false;
			sim->create_S(m_->s_);
		}
		m_->samples_list_.set_target();
	}
	if(sim->is_created()){
		sim->set_observable(which);
		sim->run(tmp_test?10:1e6,m_->tmax_);
#pragma omp critical(samples_list_)
		{
			if(m_->samples_list_.find_sorted(sim,MCSim::sort_by_param_for_merge)){
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
		std::cerr<<__PRETTY_FUNCTION__<<" : not valid parameter : "<<param<<std::endl;
		return NULL;
	}
}
/*}*/
/*}*/

/*{Minimization*/
void VMCMinimization::Minimization::set(Parseur& P, std::string& path, std::string& basename){
	unsigned int i(0);
	if(P.find("tmax",i,false)){ tmax_ = P.get<unsigned int>("tmax"); }
	else {
		tmax_ = 1;
		std::string msg("assume tmax=1s");
		std::cout<<"#"+msg<<std::endl;
		pso_info_.item(msg);
	}

	if(P.find("load",i,false)){
		Time chrono;
		std::string filename(P.get<std::string>(i));
		std::string msg("loading samples from "+filename);
		std::cout<<"#"+msg<<std::flush;

		IOFiles in(filename,false);
		std::string n_samples(load(in,path,basename));

		pso_info_.title("Minimization",'>');
		std::string msg_end(" ("+n_samples+" samples loaded in "+my::tostring(chrono.elapsed())+"s)");
		pso_info_.item(msg+msg_end);
		std::cout<<msg_end<<std::endl;
	} else {
		create(P,path,basename);
	}
}

void VMCMinimization::Minimization::create(Parseur& P, std::string& path, std::string& basename){
	tmax_= P.get<unsigned int>("tmax");
	dof_ = P.get<unsigned int>("dof");
	s_ = new System(P);
	ps_= new Vector<double>[dof_];

	/*!the next block is required to compute the J setup*/
	Vector<double> tmp;
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.set_bonds(s_);

	std::string msg("no samples loaded");
	std::cout<<"#"+msg<<std::endl;
	pso_info_.title("Minimization",'>');
	pso_info_.item(msg);

	set_phase_space(P);

	Linux command;
	path = cs.get_path();
	path+= my::tostring(dof_)+"dof/";
	basename = "-" + cs.get_filename();
	command.mkdir(path);
}

std::string VMCMinimization::Minimization::load(IOFiles& in, std::string& path, std::string& basename){
	s_ = new System(in);
	dof_ = in.read<unsigned int>();
	ps_= new Vector<double>[dof_];

	/*!the next block is required to compute the J setup*/
	Vector<double> tmp;
	CreateSystem cs(s_);
	cs.init(&tmp,NULL);
	cs.set_bonds(s_);

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
	pso_info_.text(header);
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
			pso_info_.item(msg);
			pso_info_.nl();
			pso_info_.lineblock(PS);
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
	out.add_header()->text(pso_info_.get());

	samples_list_.set_target();
	while(samples_list_.target_next()){ samples_list_.get().write(out); }
}
/*}*/
