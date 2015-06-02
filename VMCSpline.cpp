#include "VMCSpline.hpp"

VMCSpline::VMCSpline(Minimization& m):
	VMCMinimization(m,"PSpline"),
	pspline_(2)
{}

void VMCSpline::init(){
	set_time();
	m_.pso_info_.title("New PSpline run",'-');
	all_min_idx_.clear();
	pspline_.set();
	std::cout<<"#######################"<<std::endl;
	std::cout<<"#new VMCSpline"<<std::endl;
	std::cout<<"#"<<get_filename()<<std::endl;
	std::cout<<"#contains "<<m_.all_results_.size()<<" samples"<<std::endl;
}

void VMCSpline::set_ps(unsigned int const& i, Vector<double> const& ps){
	m_.set_ps(i,ps);
}

void VMCSpline::run(unsigned int const& explore_around_minima){
	if(m_.all_results_.size()){
		std::cout<<"#run"<<std::endl;
		search_minima();

		if(explore_around_minima){
#pragma omp parallel for
			for(unsigned int i=0;i<all_min_idx_.size();i++){
				std::vector<unsigned int>* min_idx(new std::vector<unsigned int>[m_.Nfreedom_]);
				for(unsigned int j(0);j<m_.Nfreedom_;j++){
					for(unsigned int k(0);k<2*explore_around_minima+1;k++){
						if(all_min_idx_[i](j)+k >= 2 && all_min_idx_[i](j)+k < m_.ps_[j].size()+2){ 
							min_idx[j].push_back(all_min_idx_[i](j)+k-2);
						}
					}
				}

				Vector<double>* local_x(new Vector<double>[m_.Nfreedom_]);
				for(unsigned int j(0);j<m_.Nfreedom_;j++){
					local_x[j].set(min_idx[j].size(),0);
					for(unsigned int k(0);k<min_idx[j].size();k++){
						local_x[j](k) = m_.ps_[j](min_idx[j][k]);
					}
				}

				Vector<unsigned int> idx(m_.Nfreedom_,0);
				while(go_through_parameter_space(local_x,idx,0,0,&VMCSpline::call_compute_vmc));
				delete[] local_x;
				delete[] min_idx;
			}
		} else {
#pragma omp parallel for
			for(unsigned int i=0;i<all_min_idx_.size();i++){
				Vector<double> param(m_.Nfreedom_);
				for(unsigned int j(0); j<m_.Nfreedom_;j++){ param(j) = m_.ps_[j](all_min_idx_[i](j)); }
				compute_vmc(param);
			}
		} 
	}
}

bool VMCSpline::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&)){
	(this->*f)(x,idx);

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<m_.Nfreedom_;i++){
			if(idx(i) == x[i].size()){ 
				if(i+1 == m_.Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	} else {
		if(idx(0) == max0){ 
			if(1 == m_.Nfreedom_){ return false; }
			idx(0) = min0;
			idx(1)++;
		}
		for(unsigned int i(1);i<m_.Nfreedom_;i++){
			if(idx(i) == x[i].size()){ 
				if(i+1 == m_.Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return go_through_parameter_space(x,idx,min0,max0,f);
}

void VMCSpline::call_compute_vmc(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_.Nfreedom_);
	for(unsigned int i(0); i<m_.Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	compute_vmc(param);
}

void VMCSpline::select_if_min(Vector<double>* x, Vector<unsigned int> const& idx){
	double f;
	double f_tmp;
	bool is_min(true);
	Vector<double> param(m_.Nfreedom_);
	for(unsigned int i(0); i<m_.Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	f = pspline_.extrapolate(param);
	if(!isnan(f)){
		for(unsigned int j(0);j<m_.Nfreedom_;j++){
			if(is_min){
				if(idx(j)+1<m_.ps_[j].size()){
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
			all_min_idx_.push_back(idx);
		}
	}
}

void VMCSpline::save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_.Nfreedom_);
	for(unsigned int i(0); i<m_.Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	(*out_)<<param<<" "<<pspline_.extrapolate(param)<<IOFiles::endl;
}

void VMCSpline::plot(){
	if(m_.Nfreedom_<4){
		out_ = new IOFiles(get_filename()+".dat",true);
		m_.all_results_.set_target();
		double min(0);
		while(m_.all_results_.target_next()){
			if(m_.all_results_.get().get_S()->get_energy().get_x()<min){ min = m_.all_results_.get().get_S()->get_energy().get_x(); }
			(*out_)<<m_.all_results_.get().get_param()<<" "<<m_.all_results_.get().get_S()->get_energy().get_x()<<IOFiles::endl;
		}
		delete out_;

		out_ = new IOFiles(get_filename()+"-spline.dat",true);
		Vector<unsigned int> idx(m_.Nfreedom_,0);
		while(go_through_parameter_space(m_.ps_,idx,0,0,&VMCSpline::save_spline_data));
		delete out_;
		out_ = NULL;

		Gnuplot plot("./",get_filename());
		plot.range("x","-2","2");
		plot.range("y","-2","2");
		plot.range("z","-1","");
		plot+="set ticslevel 0";
		if(m_.Nfreedom_==2){
			plot+="splot '"+get_filename()+"-spline.dat' u 1:2:3 lt 7 lc 1 t 'spline',\\";
			plot+="      '"+get_filename()+".dat'        u 1:2:3 lt 7 lc 2 notitle";
		} 
		if(m_.Nfreedom_==3){
			plot+="splot '"+get_filename()+"-spline.dat' u 1:2:3:($4>-0.93?0:exp("+my::tostring(min)+"-$4)) w p lt 7 lc 1 ps variable t 'spline',\\";
			plot+="      '"+get_filename()+".dat'        u 1:2:3:($4>-0.93?0:exp("+my::tostring(min)+"-$4)) w p lt 7 lc 2 ps variable notitle";
		} 
		plot.save_file();
	} else {
		std::cerr<<"void VMCSpline::plot() : how to plot "<<m_.Nfreedom_<<"-dimensional data ?"<<std::endl;
	}
}

void VMCSpline::search_minima(){
	all_min_idx_.clear();
	while(m_.all_results_.target_next()){
		pspline_.add_data(m_.all_results_.get().get_param(),m_.all_results_.get().get_S()->get_energy().get_x());
	}
	pspline_.compute_weights();

#pragma omp parallel 
	{
		unsigned int min0(omp_get_thread_num()*m_.ps_[0].size()/omp_get_num_threads());
		unsigned int max0((omp_get_thread_num()+1)*m_.ps_[0].size()/omp_get_num_threads());
		Vector<unsigned int> idx(m_.Nfreedom_,0);
		idx(0) = min0;
		while(go_through_parameter_space(m_.ps_,idx,min0,max0,&VMCSpline::select_if_min));
	}
	for(unsigned int j(0); j<all_min_idx_.size();j++){ 
		Vector<double> param(m_.Nfreedom_);
		for(unsigned int i(0); i<m_.Nfreedom_;i++){ param(i) = m_.ps_[i](all_min_idx_[j](i)); }
	}
	std::cout<<"#Minima : found "<<all_min_idx_.size()<<std::endl;
	for(unsigned int i(0); i<all_min_idx_.size();i++){
		std::cout<<"# "<<i<<" : ";
		for(unsigned int j(0);j<m_.Nfreedom_;j++){
			std::cout<<m_.ps_[j](all_min_idx_[i](j))<<" ";
		}
		std::cout<<std::endl;
	}
}

void VMCSpline::print() const {
	Vector<double> param(m_.Nfreedom_);
	for(unsigned int i(0);i<all_min_idx_.size();i++){
		for(unsigned int j(0); j<m_.Nfreedom_;j++){ param(j) = m_.ps_[j](all_min_idx_[i](j)); }
		std::cout<<param<<" "<<pspline_.extrapolate(param)<<std::endl;
	}
}
