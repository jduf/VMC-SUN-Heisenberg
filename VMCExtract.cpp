#include "VMCExtract.hpp"

VMCExtract::VMCExtract(IOFiles& in, unsigned int const& min_sort, unsigned int const& max_sort):
	VMCMinimization(in,false,"EXT")
{
	List<MCSim> possible_samples;
	unsigned int n_samples(in.read<unsigned int>());

	std::cout<<"#######################"<<std::endl;
	std::string msg("VMCExtract");
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.title(msg,'-');
	msg = "sorts all samples and discard information for the ones far from the minimum";
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.text(msg);

	std::cout<<"#"<<in.get_filename()<<std::endl;
	m_->info_.item(in.get_filename());

	unsigned int iter(0);
	if(min_sort && max_sort){
		std::shared_ptr<MCSim> sim;
		List<MCSim> sorted_samples;

		/*update sorted_samples so it always contains the min_sort lowest
		 *samples, the other ones are set in possible_samples list*/
		double Emax(0.0);
		while(iter++<n_samples){
			sim = std::make_shared<MCSim>(in);
			if(sorted_samples.size() < min_sort){
				sorted_samples.add_sort(sim,MCSim::sort_by_E);
			} else if(sim->get_energy().get_x()<Emax){
				sorted_samples.add_sort(sim,MCSim::sort_by_E);
				possible_samples.add_sort(sorted_samples.last_ptr(),MCSim::sort_for_merge);
				sorted_samples.pop_end();
				Emax = sorted_samples.last().get_energy().get_x();
			} else { possible_samples.add_end(sim); }

			if(!(iter%1000)){ my::display_progress(iter,n_samples,"loading (A): "); }
		}

		iter=0;
		/*add to m_->samples all samples of possible_samples that are close in energy*/
		double Emin(sorted_samples.first().get_energy().get_x());
		double Etmp(0.0);
		double dEtmp(0.0);
		possible_samples.set_target();
		while(possible_samples.target_next()){
			if(!(iter++%1000)){ my::display_progress(iter,possible_samples.size(),"loading (B): "); }
			Etmp = possible_samples.get().get_energy().get_x();
			dEtmp= possible_samples.get().get_energy().get_dx();
			if( Emax > Etmp-dEtmp || Emin > Etmp-2*dEtmp ){
				if(sorted_samples.size() >= max_sort){
					if( Etmp < sorted_samples.last().get_energy().get_x() ){
						possible_samples.add_sort(sorted_samples.last_ptr(),MCSim::sort_for_merge);
						sorted_samples.pop_end();

						sorted_samples.add_sort(possible_samples.get_ptr(),MCSim::sort_by_E);
						possible_samples.pop_target();
					}
				} else {
					sorted_samples.add_sort(possible_samples.get_ptr(),MCSim::sort_by_E);
					possible_samples.pop_target();
				}
			}
		}

		/*set m_->samples_ so it contains the samples of sorted_samples*/
		sorted_samples.set_target();
		while(sorted_samples.target_next()){
			m_->samples_.add_sort(sorted_samples.get_ptr(),MCSim::sort_for_merge);
		}
	} else {
		while(iter++<n_samples){ m_->samples_.add_end(std::make_shared<MCSim>(in)); }
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

	msg="keep only "+my::tostring(m_->samples_.size())+" MCSim";
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);;

	msg="and keep  "+my::tostring(dis_sim_.size())+" DiscardedSim";
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);

	dis_sim_.set_target();
	if(dis_sim_.target_next()){
		Vector<double> param(dis_sim_.get().get_param());
		while(dis_sim_.target_next()){
			if(MCSim::sort_by_param_for_merge(param,dis_sim_.get().get_param())!=1){
				std::cerr<<__PRETTY_FUNCTION__<<" : the discarded samples are not sorted correctly after construction"<<std::endl;
			}
			param = dis_sim_.get().get_param();
		}
	}

	m_->samples_.set_target();
	if(m_->samples_.target_next()){
		Vector<double> param(m_->samples_.get().get_param());
		while(m_->samples_.target_next()){
			if(MCSim::sort_by_param_for_merge(param,m_->samples_.get().get_param())!=1){
				std::cerr<<__PRETTY_FUNCTION__<<" : the samples are not sorted correctly after construction"<<std::endl;
			}
			param = m_->samples_.get().get_param();
		}
	}
}

void VMCExtract::refine(Vector<unsigned int> const& which_obs, double const& dEoE, unsigned int const& t, unsigned int maxiter){
	total_eval_ = m_->samples_.size();
	if(maxiter){ m_->tmax_ = t; }
	else {
		maxiter = std::min((unsigned int)sqrt(t/(5*total_eval_)),(unsigned int)5);
		m_->tmax_ = std::min(maxiter*5,(unsigned int)60);
	}

	if(total_eval_ && m_->tmax_){
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
		std::cerr<<__PRETTY_FUNCTION__<<" : not enough time to refine "<<my::tostring(m_->samples_.size())<<" samples (would need at least "<<my::tostring(total_eval_*5)<<"s)"<<std::endl;
		m_->samples_.set();
	}
}

void VMCExtract::save(std::string dirname) const {
	if(m_->samples_.size()){
		my::ensure_trailing_slash(dirname);
		dirname += get_path();
		set_time();
		Linux().mkpath(dirname.c_str());
		IOFiles out(dirname+get_filename()+".jdbin",true,false);
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

List<MCSim>::Node* VMCExtract::analyse(std::string const& path, std::string const& filename, List<MCSim>& kept_samples) const {
	List<MCSim>::Node* target(NULL);
	if(m_->samples_.size()){
		double E(666);
		double Erange(-666);
		double norm(0.0);
		double v_norm(0.0);
		double m_norm(0.0);
		double tmp(0.0);
		unsigned int i(0);
		Vector<double> param;
		std::set<unsigned int> n_zone;

		m_->samples_.set_target();
		m_->samples_.target_next();
		param = m_->samples_.get().get_param();
		IOFiles data_n(path+filename+"-n.dat",true,false);
		/*!compute the mean distance between paramters set (norm) and related
		 * variance*/
		do {
			norm = sqrt((param-m_->samples_.get().get_param()).norm_squared())/param.size();
			param = m_->samples_.get().get_param();
			data_n<<i<<" "<<norm<<IOFiles::endl;

			i++;
			tmp = norm-m_norm;
			m_norm += tmp/i;
			v_norm += tmp*(norm-m_norm);

			if(E>m_->samples_.get().get_energy().get_x()){
				E = m_->samples_.get().get_energy().get_x();
				target = m_->samples_.get_target();
			}
		} while(m_->samples_.target_next());
		v_norm = sqrt(v_norm/(i-1));

		/*!split the whole list of samples into intervals of size proportional
		 * to the variance*/
		unsigned int k(0);
		while(n_zone.size()<50 && k++<5){
			m_->samples_.set_target();
			m_->samples_.target_next();
			param = m_->samples_.get().get_param();
			E = m_->samples_.get().get_energy().get_x();
			Erange = E;
			i=0;
			do {
				Erange = std::max(Erange,m_->samples_.get().get_energy().get_x()+3*m_->samples_.get().get_energy().get_dx());
				norm = sqrt((param-m_->samples_.get().get_param()).norm_squared())/param.size();
				param = m_->samples_.get().get_param();

				if(norm>m_norm+2*v_norm/k){ n_zone.insert(i); }
				i++;
			} while(m_->samples_.target_next());
		}
		for(unsigned int i(0);i<5;i++){ n_zone.insert((i+1)*m_->samples_.size()/5); }

		/*!in each interval, find and save the local minima*/
		if(n_zone.size()){
			i=0;
			unsigned int idx(0);
			unsigned int j(0);
			m_->samples_.set_target();
			std::vector<unsigned int> n;
			std::set<unsigned int>::iterator z(n_zone.begin());
			while(m_->samples_.target_next()){
				if(E>m_->samples_.get().get_energy().get_x()){
					E = m_->samples_.get().get_energy().get_x();
					idx = i;
				}
				if(++i==*z){
					z++;
					n.push_back(idx);
					E=0.0;
				}
			}

			dis_sim_.set_target();
			dis_sim_.target_next();
			m_->samples_.set_target();
			m_->samples_.target_next();
			param = dis_sim_.get().get_param();
			IOFiles data_Er(path+filename+"-Er.dat",true,false);
			data_Er.precision(15);
			bool remain_samples(true);
			bool keepon;

			j=0;
			i=0;
			norm=0.0;
			do{
				if(remain_samples && MCSim::sort_by_param_for_merge(m_->samples_.get().get_param(),dis_sim_.get().get_param())){
					norm += (param-m_->samples_.get().get_param()).norm_squared();
					if(i++!=n[j]){
						data_Er<<norm<<" "<<m_->samples_.get().get_energy()<<" 1"<<IOFiles::endl;
					} else {
						data_Er<<norm<<" "<<m_->samples_.get().get_energy()<<" 2"<<IOFiles::endl;
						kept_samples.add_sort(m_->samples_.get_ptr(),MCSim::sort_by_E);
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
			} while(keepon);

			Gnuplot gp(path,filename);
			gp.range("y","",Erange);
			gp+="plot '"+filename+"-Er.dat' u ($6==0?$1:1/0):2:3 w e notitle,\\";
			gp+="     '"+filename+"-Er.dat' u ($6==1?$1:1/0):2:3 w e notitle,\\";
			gp+="     '"+filename+"-Er.dat' u ($6==2?$1:1/0):2:3 w e notitle";
			gp.save_file();
			gp.create_image(true,"png");
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : can't find characteristic spacing"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples"<<std::endl; }
	return target;
}

unsigned int VMCExtract::DiscardedSim::sort(DiscardedSim const& list, DiscardedSim const& new_elem){
	return MCSim::sort_by_param_for_merge(list.param_, new_elem.param_);
}
