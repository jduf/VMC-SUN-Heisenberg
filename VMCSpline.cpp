#include "VMCSpline.hpp"

VMCSpline::VMCSpline(VMCMinimization const& m):
	VMCMinimization(m,"PSpline"),
	pspline_(2)
{}

/*{Public methods*/
void VMCSpline::init(){
	set_time();
	pspline_.set();
	list_min_idx_.clear();

	std::cout<<"#######################"<<std::endl;
	std::cout<<"#initialize VMCSpline"<<std::endl;
	std::cout<<"#"<<get_filename()<<std::endl;
	m_->pso_info_.title("New PSpline run",'-');
	m_->pso_info_.item(get_filename());

	if(m_->samples_list_.size()){ search_minima(); }
	else { 
		std::cerr<<"VMCSpline::run() : empty samples_list_"<<std::endl; 
		m_->pso_info_.item("error : empty samples_list_");
	}
}

void VMCSpline::run(unsigned int const& explore_around_minima){
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
				while(go_through_parameter_space(local_x,idx,0,0,&VMCSpline::evaluate));
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

void VMCSpline::plot(){
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
		while(go_through_parameter_space(m_->ps_,idx,0,0,&VMCSpline::save_spline_data));
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
		std::cerr<<"void VMCSpline::plot() : how to plot "<<m_->Nfreedom_<<"-dimensional data ?"<<std::endl;
	}
}

void VMCSpline::print(){
	Vector<double> param(m_->Nfreedom_);
#pragma omp parallel for firstprivate(param)
	for(unsigned int i=0;i<list_min_idx_.size();i++){
		for(unsigned int j(0);j<m_->Nfreedom_;j++){ param(j) = m_->ps_[j](list_min_idx_[i](j)); }
		std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(param));
		if(sim.get()){
#pragma omp critical
			std::cout<<param<<" "<<pspline_.extrapolate(param)<<" "<<sim->get_S()->get_energy()<<std::endl;
		}
	}
}
/*}*/

/*{Private methods*/
void VMCSpline::search_minima(){
	list_min_idx_.clear();
	pspline_.set(m_->samples_list_.size());
	unsigned int i(0);
	while(m_->samples_list_.target_next()){
		pspline_.add_data(i++,m_->samples_list_.get().get_param(),m_->samples_list_.get().get_S()->get_energy().get_x());
	}
	Time chrono;
	std::string msg1("computing the weights");
	std::cout<<"#"<<msg1<<std::flush;

	pspline_.compute_weights();

	std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s)");
	std::cout<<msg2<<std::endl;
	m_->pso_info_.item(msg1+msg2);
	chrono.set();
	msg1="searching for minima";
	std::cout<<"#"<<msg1<<std::flush;

#pragma omp parallel 
	{
		unsigned int min0(omp_get_thread_num()*m_->ps_[0].size()/omp_get_num_threads());
		unsigned int max0((omp_get_thread_num()+1)*m_->ps_[0].size()/omp_get_num_threads());
		Vector<unsigned int> idx(m_->Nfreedom_,0);
		idx(0) = min0;
		while(go_through_parameter_space(m_->ps_,idx,min0,max0,&VMCSpline::select_if_min));
	}
	for(unsigned int j(0); j<list_min_idx_.size();j++){ 
		Vector<double> param(m_->Nfreedom_);
		for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = m_->ps_[i](list_min_idx_[j](i)); }
	}

	msg2=" (found "+my::tostring(list_min_idx_.size())+" minimas in "+my::tostring(chrono.elapsed())+"s)";
	std::cout<<msg2<<std::endl;
	m_->pso_info_.item(msg1+msg2);
}

bool VMCSpline::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&)){
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

void VMCSpline::select_if_min(Vector<double>* x, Vector<unsigned int> const& idx){
	double f;
	double f_tmp;
	bool is_min(true);
	Vector<double> param(m_->Nfreedom_);
	for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	f = pspline_.extrapolate(param);
	if(!isnan(f)){
		for(unsigned int j(0);j<m_->Nfreedom_;j++){
			if(is_min){
				if(idx(j)+1<m_->ps_[j].size()){
					param(j) = x[j](idx(j)+1);
					f_tmp = pspline_.extrapolate(param);
					if(isnan(f_tmp) || f>f_tmp){ is_min = false; }
				}
				if(is_min && idx(j)>0){
					param(j) = x[j](idx(j)-1);
					f_tmp = pspline_.extrapolate(param);
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

void VMCSpline::save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->Nfreedom_);
	for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	(*out_)<<param<<" "<<pspline_.extrapolate(param)<<IOFiles::endl;
}

void VMCSpline::evaluate(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->Nfreedom_);
	for(unsigned int i(0); i<m_->Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	VMCMinimization::evaluate(param);
}
/*}*/
