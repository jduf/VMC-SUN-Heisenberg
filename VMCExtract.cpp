#include "VMCExtract.hpp"

VMCExtract::VMCExtract(IOFiles& in):
	VMCMinimization(in,false,"EXT")
{
	unsigned int n_samples(in.read<unsigned int>());
	unsigned int iter(0);
	std::cout<<"#######################"<<std::endl;
	std::string msg("VMCExtract");
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.title(msg,'-');

	std::cout<<"#"<<in.get_filename()<<std::endl;
	m_->info_.item(in.get_filename());

	std::shared_ptr<MCSim> sim;
	List<MCSim> sorted_samples;
	while(iter++<100 && iter<n_samples){
		sim = std::make_shared<MCSim>(in);
		sorted_samples.add_sort(sim,MCSim::sort_by_E);
	}

	List<MCSim> possible_samples;
	double Emax(sorted_samples.last().get_energy().get_x());
	Vector<double> param(sim->get_param());
	iter--; //otherwise one sample is skipped
	while(iter++<n_samples){
		sim = std::make_shared<MCSim>(in);

		if(MCSim::sort_by_param_for_merge(param,sim->get_param())==1){
			if(sim->get_energy().get_x()<Emax){
				sorted_samples.add_sort(sim,MCSim::sort_by_E);
				possible_samples.add_sort(sorted_samples.last_ptr(),MCSim::sort_for_merge);
				sorted_samples.pop_end();
				Emax = sorted_samples.last().get_energy().get_x();
			} else { possible_samples.add_end(sim); }

			param = sim->get_param();
		} else {
			iter = n_samples;
			std::cerr<<__PRETTY_FUNCTION__<<" : the samples are not sorted correctly"<<std::endl;
		}
		if(!(iter%1000)){ my::display_progress(iter,n_samples,"loading : "); }
	}

	sorted_samples.set_target();
	while(sorted_samples.target_next()){
		m_->samples_.add_sort(sorted_samples.get_ptr(),MCSim::sort_for_merge);
	}
	
	double Emin(sorted_samples.first().get_energy().get_x());
	possible_samples.set_target();
	while(possible_samples.target_next()){
		if(
				Emax > possible_samples.get().get_energy().get_x()-possible_samples.get().get_energy().get_dx()
				||
				Emin > possible_samples.get().get_energy().get_x()-2*possible_samples.get().get_energy().get_dx()
		  )
		{
			m_->samples_.add_sort(possible_samples.get_ptr(),MCSim::sort_for_merge);
			possible_samples.pop_target();
		}
	}

	load_filenames(in);

	msg="contains "+my::tostring(n_samples)+" samples";
	if(in.get_filename().find("EXT") != std::string::npos){
		in>>n_samples;
		for(unsigned int i(0);i<n_samples;i++){
			dis_sim_.add_end(std::make_shared<VMCExtract::DiscardedSim>(in));
			if(!(i%1000)){ my::display_progress(i,n_samples,"loading : "); }
		}

		possible_samples.set_target();
		while(possible_samples.target_next()){
			dis_sim_.add_sort(std::make_shared<DiscardedSim>(possible_samples.get_ptr()),VMCExtract::DiscardedSim::sort);
		}

		msg+=" and "+my::tostring(n_samples)+" discarded samples";
	} else {
		possible_samples.set_target();
		while(possible_samples.target_next()){
			dis_sim_.add_end(std::make_shared<DiscardedSim>(possible_samples.get_ptr()));
		}
	}
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);

	msg="keep only "+my::tostring(m_->samples_.size())+" samples";
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);

	dis_sim_.set_target();
	dis_sim_.target_next();
	param = dis_sim_.get().get_param();
	while(dis_sim_.target_next()){
		if(!MCSim::sort_by_param_for_merge(param,dis_sim_.get().get_param())==1){
			std::cerr<<__PRETTY_FUNCTION__<<" : the samples are not sorted correctly after construction"<<std::endl;
		}
		param = dis_sim_.get().get_param();
	}

	m_->samples_.set_target();
	m_->samples_.target_next();
	param = m_->samples_.get().get_param();
	while(m_->samples_.target_next()){
		if(!MCSim::sort_by_param_for_merge(param,m_->samples_.get().get_param())==1){
			std::cerr<<__PRETTY_FUNCTION__<<" : the samples are not sorted correctly after construction"<<std::endl;
		}
		param = m_->samples_.get().get_param();
	}
}

void VMCExtract::refine(Vector<unsigned int> const& which_obs, double const& dEoE, unsigned int const& ttotal){
	unsigned int maxiter(sqrt(ttotal/5/m_->samples_.size()));
	m_->tmax_ = std::min(maxiter*5,(unsigned int)60);
	maxiter = std::min(maxiter,(unsigned int)5);

	if(m_->samples_.size() && m_->tmax_){
		total_eval_ = m_->samples_.size();
		progress_ = 0;

		std::string msg("refines "+my::tostring(total_eval_)+" samples (max time "+my::tostring(total_eval_*m_->tmax_*maxiter)+"s)");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = RST::math("t_{max} = "+my::tostring(m_->tmax_)+"s")+", "+RST::math("\\mathrm{d}E/E="+my::tostring(dEoE)) + " , maxiter="+my::tostring(maxiter);
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		m_->samples_.set_target();
		while(m_->samples_.target_next()){ evaluate_until_precision(m_->samples_.get().get_param(),which_obs,dEoE,maxiter); }
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : not enough time to refine "<<my::tostring(m_->samples_.size())<<" samples (would need "<<my::tostring(total_eval_*m_->tmax_*maxiter)<<")"<<std::endl;
		m_->samples_.set();
	}
}

void VMCExtract::save(std::string const& path) const {
	if(m_->samples_.size()){
		set_time();
		Linux().mkpath((path+get_path()).c_str());
		IOFiles out(path+get_path()+get_filename()+".jdbin",true);
		out.add_header()->text(RST::textbf("Contains the "+my::tostring(m_->samples_.size())+" best samples and a list "+my::tostring(dis_sim_.size())+ " of discarded samples"));
		VMCMinimization::save(out);
		out<<dis_sim_.size();

		dis_sim_.set_target();
		while(dis_sim_.target_next()){ dis_sim_.get().write(out); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no sample to save"<<std::endl; }
}

void VMCExtract::print() const {
	m_->samples_.set_target();
	while(m_->samples_.target_next()){
		std::cout<<m_->samples_.get().get_param()<<std::endl;
		m_->samples_.get().print(0);
	}
}

void VMCExtract::plot() const {
	std::string path("./");
	std::string filename("test");
	if(m_->samples_.size()){
		double E;
		double Erange;
		double norm(0.0);
		unsigned int i(0);
		unsigned int j(0);
		Vector<double> param;
		std::vector<unsigned int> n;

		m_->samples_.set_target();
		m_->samples_.target_next();
		param = m_->samples_.get().get_param();
		E = m_->samples_.get().get_energy().get_x();
		Erange = E;
		IOFiles data_n(path+filename+"-n.dat",true);
		do {
			Erange = std::max(Erange,m_->samples_.get().get_energy().get_x()+10*m_->samples_.get().get_energy().get_dx());
			norm = sqrt((param-m_->samples_.get().get_param()).norm_squared())/param.size();
			param = m_->samples_.get().get_param();

			data_n<<i<<" "<<norm<<IOFiles::endl;
			if(norm>0.1){
				n.push_back(i); 
				std::cout<<i<<" | "<<param<<std::endl;
			}
			i++;
		} while(m_->samples_.target_next());

		i=0;
		unsigned int idx(0);
		m_->samples_.set_target();
		while(m_->samples_.target_next()){
			i++; 
			if(E>m_->samples_.get().get_energy().get_x()){
				E=m_->samples_.get().get_energy().get_x();
				idx = i;
			}
			if(i==n[j]){
				n[j] = idx;
				j++;
				E=0.0;
			}
		}

		dis_sim_.set_target();
		dis_sim_.target_next();
		m_->samples_.set_target();
		m_->samples_.target_next();
		param = dis_sim_.get().get_param();
		IOFiles data_Er(path+filename+"-Er.dat",true);
		data_Er.precision(15);
		bool remain_samples(true);
		bool keepon;

		j=0;
		i=0;
		norm=0.0;
		do{
			if(remain_samples && MCSim::sort_by_param_for_merge(m_->samples_.get().get_param(),dis_sim_.get().get_param())){
				norm += (param-m_->samples_.get().get_param()).norm_squared();
				if(++i!=n[j]){
					data_Er<<norm<<" "<<m_->samples_.get().get_energy()<<" 1"<<IOFiles::endl;
				} else {
					data_Er<<norm<<" "<<m_->samples_.get().get_energy()<<" 2"<<IOFiles::endl;
					j++;
				}

				param = m_->samples_.get().get_param();
				remain_samples = m_->samples_.target_next();
			} else {
				norm += (param-dis_sim_.get().get_param()).norm_squared();
				data_Er<<norm<<" "<<dis_sim_.get().get_energy()<<" 0"<<IOFiles::endl;

				param = dis_sim_.get().get_param();
				keepon = dis_sim_.target_next() ;
			}
			if(my::are_equal(norm,floor(norm),1e-15,1e-15)){ std::cerr<<"strange"<<std::endl; }
		} while(keepon);

		Gnuplot gp(path,filename);
		gp.range("y","",Erange);
		gp+="plot '"+filename+"-Er.dat' u ($6==0?$1:1/0):2:3 w e notitle,\\";
		gp+="     '"+filename+"-Er.dat' u ($6==1?$1:1/0):2:3 w e notitle,\\";
		gp+="     '"+filename+"-Er.dat' u ($6==2?$1:1/0):2:3 w e notitle";
		gp.save_file();
		gp.create_image(true,true);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples"<<std::endl; }
}

unsigned int VMCExtract::DiscardedSim::sort(DiscardedSim const& list, DiscardedSim const& new_elem){
	return MCSim::sort_by_param_for_merge(list.param_, new_elem.param_);
}
