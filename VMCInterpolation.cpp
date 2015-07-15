#include "VMCInterpolation.hpp"

VMCInterpolation::VMCInterpolation(VMCMinimization const& m):
	VMCMinimization(m,"I"),
	interp_(m_->Nfreedom_,pow(m_->ps_size_/m_->samples_list_.size(),1./m_->Nfreedom_))
{
	interp_.select_basis_function(7);
	std::cerr<<"VMCInterpolation::VMCInterpolation(VMCMinimization const& m): no 'dx' set for Interpolation"<<std::endl;
}

/*{public methods*/
void VMCInterpolation::init(){
	set_time();
	interp_.set();
	list_min_idx_.clear();

	std::cout<<"#######################"<<std::endl;
	std::string msg("new VMCInterpolation");
	std::cout<<"#"<<msg<<std::endl;
	m_->pso_info_.title(msg,'-');

	std::cout<<"#"<<get_filename()<<std::endl;
	m_->pso_info_.item(get_filename());

	msg="contains "+my::tostring(m_->samples_list_.size())+" samples";
	std::cout<<"#"<<msg<<std::endl;
	m_->pso_info_.item(msg);

	if(m_->samples_list_.size()){ search_minima(); }
	else {
		std::cerr<<"VMCInterpolation::run() : empty samples_list_"<<std::endl;
		m_->pso_info_.item("error : empty samples_list_");
	}
}

void VMCInterpolation::run(unsigned int const& explore_around_minima){
	unsigned int lmis(list_min_idx_.size());
	if(lmis){
		unsigned int nsim(pow(2*explore_around_minima+1,m_->Nfreedom_));
		std::string msg1("measuring "+my::tostring(nsim)+" samples per minima (time estimated "+my::tostring(lmis*m_->tmax_*nsim/omp_get_max_threads())+"s");
		std::cout<<"#"<<msg1<<std::flush;
		Time chrono;

		if(explore_around_minima){
#pragma omp parallel for
			for(unsigned int i=0;i<lmis;i++){
				std::vector<unsigned int>* min_idx(new std::vector<unsigned int>[m_->Nfreedom_]);
				for(unsigned int j(0);j<m_->Nfreedom_;j++){
					for(unsigned int k(0);k<2*explore_around_minima+1;k++){
						if(list_min_idx_[i](j)+k >= 2 && list_min_idx_[i](j)+k < m_->ps_[j].size()+2){
							min_idx[j].push_back(list_min_idx_[i](j)+k-2);
						}
					}
				}

				Vector<double>* local_x(new Vector<double>[m_->Nfreedom_]);
				for(unsigned int j(0);j<m_->Nfreedom_;j++){
					local_x[j].set(min_idx[j].size(),0);
					for(unsigned int k(0);k<min_idx[j].size();k++){
						local_x[j](k) = m_->ps_[j](min_idx[j][k]);
					}
				}

				Vector<unsigned int> idx(m_->Nfreedom_,0);
				while(go_through_parameter_space(local_x,idx,0,0,&VMCInterpolation::evaluate));
				delete[] local_x;
				delete[] min_idx;
			}
		} else {
#pragma omp parallel for
			for(unsigned int i=0;i<lmis;i++){
				evaluate(m_->ps_,list_min_idx_[i]);
			}
		}

		std::string msg2(", done in "+my::tostring(chrono.elapsed())+"s)");
		std::cout<<msg2<<std::endl;
		m_->pso_info_.item(msg1+msg2);
	}
}

void VMCInterpolation::plot(){
	if(m_->Nfreedom_<4){
		out_ = new IOFiles(get_filename()+".dat",true);
		m_->samples_list_.set_target();
		double min(0);
		while(m_->samples_list_.target_next()){
			if(m_->samples_list_.get().get_S()->get_energy().get_x()<min){ min = m_->samples_list_.get().get_S()->get_energy().get_x(); }
			(*out_)<<m_->samples_list_.get().get_param()<<" "<<m_->samples_list_.get().get_S()->get_energy().get_x()<<IOFiles::endl;
		}
		delete out_;

		out_ = new IOFiles(get_filename()+"-spline.dat",true);
		Vector<unsigned int> idx(m_->Nfreedom_,0);
		while(go_through_parameter_space(m_->ps_,idx,0,0,&VMCInterpolation::save_interp_data));
		delete out_;
		out_ = NULL;

		Gnuplot plot("./",get_filename());
		plot.range("x","-2","2");
		plot.range("y","-2","2");
		plot.range("z","-1","");
		plot+="set ticslevel 0";
		if(m_->Nfreedom_==2){
			plot+="splot '"+get_filename()+"-spline.dat' u 1:2:3 lt 7 lc 1 t 'spline',\\";
			plot+="      '"+get_filename()+".dat'        u 1:2:3 lt 7 lc 2 notitle";
		}
		if(m_->Nfreedom_==3){
			plot+="splot '"+get_filename()+"-spline.dat' u 1:2:3:($4>-0.93?0:exp("+my::tostring(min)+"-$4)) w p lt 7 lc 1 ps variable t 'spline',\\";
			plot+="      '"+get_filename()+".dat'        u 1:2:3:($4>-0.93?0:exp("+my::tostring(min)+"-$4)) w p lt 7 lc 2 ps variable notitle";
		}
		plot.save_file();
	} else {
		std::cerr<<"void VMCInterpolation::plot() : how to plot "<<m_->Nfreedom_<<"-dimensional data ?"<<std::endl;
	}
}

void VMCInterpolation::print(){
	Vector<double> param(m_->Nfreedom_);
#pragma omp parallel for firstprivate(param)
	for(unsigned int i=0;i<list_min_idx_.size();i++){
		for(unsigned int j(0);j<m_->Nfreedom_;j++){ param(j) = m_->ps_[j](list_min_idx_[i](j)); }
		std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(param));
		if(sim.get()){
#pragma omp critical
			std::cerr<<param<<" "<<interp_.extrapolate(param)<<" "<<sim->get_S()->get_energy()<<std::endl;
		}
	}
}
/*}*/

/*{private methods*/
void VMCInterpolation::search_minima(){
	list_min_idx_.clear();
	interp_.set(m_->samples_list_.size());
	unsigned int i(0);
	while(m_->samples_list_.target_next()){
		interp_.add_data(i++,m_->samples_list_.get().get_param(),m_->samples_list_.get().get_S()->get_energy().get_x());
	}
	Time chrono;
	std::string msg1("computing the weights");
	std::cout<<"#"<<msg1<<std::flush;

	interp_.compute_weights();

	std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s)");
	std::cout<<msg2<<std::endl;
	m_->pso_info_.item(msg1+msg2);
	chrono.set();
	msg1="search for minima";
	if(m_->ps_size_ < 1e5){
		msg1="exhaustive "+msg1;
		std::cout<<"#"<<msg1<<std::flush;
#pragma omp parallel
		{
			unsigned int min0(omp_get_thread_num()*m_->ps_[0].size()/omp_get_num_threads());
			unsigned int max0((omp_get_thread_num()+1)*m_->ps_[0].size()/omp_get_num_threads());
			Vector<unsigned int> idx(m_->Nfreedom_,0);
			idx(0) = min0;
			while(go_through_parameter_space(m_->ps_,idx,min0,max0,&VMCInterpolation::select_if_min));
		}
	} else {
		msg1="random "+msg1;
		std::cout<<"#"<<msg1<<std::flush;
		RandArray<unsigned int>* rnd(new RandArray<unsigned int>[m_->Nfreedom_]);
		unsigned int max_threads(omp_get_max_threads());
		for(unsigned int i(0);i<m_->Nfreedom_;i++){
			rnd[i].set(max_threads);
			for(unsigned int j(0);j<max_threads;j++){
				rnd[i].set(j,0,m_->ps_[i].size()-1);
			}
		}
#pragma omp parallel for
		for(unsigned int i=0;i<100;i++){
			Vector<double> param(m_->Nfreedom_);
			Vector<unsigned int> pidx(m_->Nfreedom_);
			unsigned int thread(omp_get_thread_num());
			for(unsigned int j(0);j<m_->Nfreedom_;j++){
				pidx(j) = rnd[j].get(thread);
				param(j) = m_->ps_[j](pidx(j));
			}

			double tmp(interp_.extrapolate(param));
			double min(tmp);
			unsigned int dir;
			unsigned int const max_step(1e3);
			for(unsigned int j(0);j<max_step;j++){
				dir = 2*m_->Nfreedom_;
				for(unsigned int k(0);k<m_->Nfreedom_;k++){
					if(pidx(k)+1<m_->ps_[k].size()){
						param(k) = m_->ps_[k](pidx(k)+1);
						tmp = interp_.extrapolate(param);
						if(!isnan(tmp) && tmp<min){
							min = tmp;
							dir = 2*k;
						}
					}
					if(pidx(k)>0){
						param(k) = m_->ps_[k](pidx(k)-1);
						tmp = interp_.extrapolate(param);
						if(!isnan(tmp) && tmp<min){
							min = tmp;
							dir = 2*k+1;
						}
					}
					param(k) = m_->ps_[k](pidx(k));
				}
				if(dir != 2*m_->Nfreedom_){
					pidx(dir/2) += (dir%2?-1:1);
					param(dir/2) = m_->ps_[dir/2](pidx(dir/2));
				} else {
					j=max_step; 
#pragma omp critical
					{
						bool add_to_list_min(true);
						for(unsigned int k(0);k<list_min_idx_.size();k++){
							if(my::are_equal(pidx,list_min_idx_[k])){
								add_to_list_min = false;
								k = list_min_idx_.size();
							}
						}
						if(add_to_list_min){ list_min_idx_.push_back(pidx); }
					}
				}
			}
		}
		delete[] rnd;
		rnd = NULL;
	}

	msg2=" (found "+my::tostring(list_min_idx_.size())+" minimas in "+my::tostring(chrono.elapsed())+"s)";
	std::cout<<msg2<<std::endl;
	m_->pso_info_.item(msg1+msg2);
}

bool VMCInterpolation::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCInterpolation::*f)(Vector<double>*, Vector<unsigned int> const&)){
	(this->*f)(x,idx);

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<m_->Nfreedom_;i++){
			if(idx(i) == x[i].size()){
				if(i+1 == m_->Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	} else {
		if(idx(0) == max0){
			if(1 == m_->Nfreedom_){ return false; }
			idx(0) = min0;
			idx(1)++;
		}
		for(unsigned int i(1);i<m_->Nfreedom_;i++){
			if(idx(i) == x[i].size()){
				if(i+1 == m_->Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return go_through_parameter_space(x,idx,min0,max0,f);
}

void VMCInterpolation::select_if_min(Vector<double>* x, Vector<unsigned int> const& idx){
	double f;
	double f_tmp;
	bool is_min(true);
	Vector<double> param(m_->Nfreedom_);
	for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	f = interp_.extrapolate(param);
	if(!isnan(f)){
		for(unsigned int j(0);j<m_->Nfreedom_;j++){
			if(is_min){
				if(idx(j)+1<m_->ps_[j].size()){
					param(j) = x[j](idx(j)+1);
					f_tmp = interp_.extrapolate(param);
					if(isnan(f_tmp) || f>f_tmp){ is_min = false; }
				}
				if(is_min && idx(j)>0){
					param(j) = x[j](idx(j)-1);
					f_tmp = interp_.extrapolate(param);
					if(isnan(f_tmp) || f>f_tmp){ is_min = false; }
				}
				param(j) = x[j](idx(j));
			}
		}
	} else { is_min = false; }

	if(is_min){
#pragma omp critical(list_min)
		{
			list_min_idx_.push_back(idx);
		}
	}
}

void VMCInterpolation::save_interp_data(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->Nfreedom_);
	for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	(*out_)<<param<<" "<<interp_.extrapolate(param)<<IOFiles::endl;
}

void VMCInterpolation::evaluate(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->Nfreedom_);
	for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(param));
	if(sim.get()){
		sim->get_S()->get_energy().complete_analysis(1e-5);
		if(std::abs(interp_.extrapolate(param)-sim->get_S()->get_energy().get_x())<sim->get_S()->get_energy().get_dx()){
#pragma omp critical
			std::cerr<<param<<" : "<<interp_.extrapolate(param)<<" <-> "<<sim->get_S()->get_energy()<<std::endl;
		}
	}
}
/*}*/
