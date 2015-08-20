#include "VMCMinimization.hpp"

/*{VMCMinimization*/
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
		std::cerr<<"VMCMinimization::VMCMinimization(Parseur const& P) : something went wrong"<<std::endl;
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
				if(m_->samples_list_.find_sorted(sim,MCSim::cmp_for_merge)){ 
					m_->samples_list_.merge_with_target(sim,MCSim::merge); 
					m_->samples_list_.get().get_S()->get_energy().complete_analysis(1e-5);
				} else {
					std::cerr<<"void VMCMinimization::refine(double const& E, double const& dE) : can't find sim in m_->samples_list_"<<std::endl;
				}
			}
		}
	} else {
		std::cerr<<"void VMCMinimization::refine(unsigned int const& Nrefine, double const& converged_criterion, unsigned int const& tmax) : there is no data"<<std::endl;
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

void VMCMinimization::save_best(unsigned int const& nsave, IOFiles& w) const {
	if(m_->samples_list_.size()){
		List<MCSim> best;
		while(m_->samples_list_.target_next()){ best.add_sort(m_->samples_list_.get_ptr(),MCSim::compare); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(w); }
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

void VMCMinimization::plot(std::string path, std::string filename) const {
	if(path==""){ path = path_; }
	if(filename==""){ filename =  get_filename(); }

	IOFiles data(path+filename+".dat",true);

	double E(0);
	double tmp;
	Vector<double> param;

	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		tmp = m_->samples_list_.get().get_S()->get_energy().get_x();
		if(tmp<E){
			E=tmp;
			param = m_->samples_list_.get().get_param();
		}
	}

	m_->samples_list_.set_target();
	List<std::pair<double,double> > r;
	std::shared_ptr<std::pair<double,double> > r_tmp;

	auto sort_by_r = [](std::pair<double,double> const& a, std::pair<double,double> const& b){
		if(my::are_equal(a.first,b.first)){ return 2; }
		if(a.first>b.first){ return 0; }
		if(a.first<b.first){ return 1; }
		return 2;
	};
	auto replace_E = [](std::pair<double,double>& a, std::pair<double,double>& b){ a.second = b.second; };

	while(m_->samples_list_.target_next()){
		r_tmp=std::make_shared<std::pair<double, double> >((param-m_->samples_list_.get().get_param()).norm_squared(),m_->samples_list_.get().get_S()->get_energy().get_x());
		data<<m_->samples_list_.get().get_param()<<" "<<r_tmp->first<<" "<<m_->samples_list_.get().get_S()->get_energy()<<IOFiles::endl;
		if(r.find_sorted(r_tmp,sort_by_r)){
			if(r.get().second > r_tmp->second){ r.merge_with_target(r_tmp,replace_E); }
		} else { r.add_after_target(r_tmp); }
	}

	List<std::pair<double,double> > r_cpy;
	bool keep_r;
	unsigned int ao(50);
	Interpolation interp_Er(1);
	interp_Er.select_basis_function(7);
	for(unsigned int i(0);i<r.size();i++){
		tmp = r[i].first;
		E = r.get().second;
		keep_r = true;
		for(unsigned int j(i>ao?i-ao:1);j<i+ao && j<r.size();j++){
			if(r[j].second<E){
				keep_r = false; 
				j = r.size();
			}
		}
		if(keep_r){
			r_cpy.add_end(std::make_shared<std::pair<double, double> >(tmp,E)); 
			interp_Er.add_data(tmp,E);
		}
	}
	double dx(1);
	interp_Er.compute_weights(dx,r_cpy.size());
	std::cout<<"ok "<<dx<<std::endl;

	IOFiles data_Er(path+filename+"-Er.dat",true);
	r_cpy.set_target();
	while(r_cpy.target_next()){
		data_Er<<r_cpy.get().first<<" "<<r_cpy.get().second<<IOFiles::endl;
	}

	Vector<double> range(0,18,0.05);
	Vector<double> Er(range.size());
	for(unsigned int i(0);i<range.size();i++){
		Er(i) = interp_Er(Vector<double>(1,range(i)));
	}

	for(unsigned int i(1);i<Er.size()-1;i++){
		if( Er(i-1)>Er(i) && Er(i)<Er(i+1) && Er(i) < -0.693) { std::cout<<range(i)<<std::endl; }
	}

	IOFiles data_interp(path+filename+"-interp.dat",true);
	for(unsigned int i(0);i<range.size();i++){
		data_interp<<range(i)<<" "<<Er(i)<<IOFiles::endl;
	}

	m_->samples_list_.target_next();
	unsigned int N(m_->samples_list_.get().get_param().size());

	Gnuplot gp(path,filename);
	gp+="Em="+my::tostring(E*0.99);
	gp.multiplot();
	gp.range("x","[:Em] writeback");
	gp.margin("0.1","0.9","0.5","0.10");
	gp+="plot '"+filename+".dat' u "+my::tostring(N+2)+":"+my::tostring(N+1)+":"+my::tostring(N+3)+" w xe notitle,\\";
	gp+="     '"+filename+"-Er.dat' u 2:1 notitle,\\";
	gp+="     '"+filename+"-interp.dat' u 2:1 w l  notitle";
	gp.margin("0.1","0.9","0.9","0.50");
	gp.tics("x");
	gp.range("x","restore");
	gp.key("left Left");
	for(unsigned int i(0);i<N;i++){
		gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(N+2)+":"+my::tostring(i+1)+":"+my::tostring(N+3)+" w xe t '$"+my::tostring(i)+"$'"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true,true);
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
						std::cerr<<"void VMCMinimization::Minimization::set_phase_space(Parseur const& P) : each range must be given this way [min:dx:max]"<<std::endl;
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
			std::cerr<<"void VMCMinimization::Minimization::set_phase_space(Parseur const& P) : provide "<<dof_<<" ranges and remove any blank space and EOL at the EOF"<<std::endl;
		}
	} else {
		std::cerr<<"void VMCMinimization::Minimization::set_phase_space(Parseur const& P) : need to provide a file containing the phase space"<<std::endl;
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
