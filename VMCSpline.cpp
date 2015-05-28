#include "VMCSpline.hpp"

VMCSpline::VMCSpline(Parseur& P):
	VMCMinimization(P,"PSpline"),
	pspline_(2)
{}

void VMCSpline::init(bool border){
	set_time();
	pso_info_.title("New PSpline run",'-');
	if(border){
		pso_info_.text("set border"+RST::nl_);
	}

#pragma omp parallel
	{
		unsigned int min0(omp_get_thread_num()*x_[0].size()/8);
		unsigned int max0((omp_get_thread_num()+1)*x_[0].size()/8);
		Vector<unsigned int> idx;
		idx.set(Nfreedom_,0);
		idx(0) = min0;
		while(set_parameter_space(idx,min0,max0,border));
	}
	unsigned int incr(border_.size()/12);
#pragma omp parallel for
	for(unsigned int i=0;i<border_.size();i+=incr){
		std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(border_[i]));
		sim->create_S(&system_param_);
		if(sim->is_created()){
			sim->run(1e6,tmax_);
#pragma omp critical(all_results_)
			{
				if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){
					all_results_.fuse_with_target(sim,MCSim::fuse);
				} else {
					all_results_.add_after_target(sim);
				}
			}
		}
	}
}

void VMCSpline::run(){
	if(all_results_.size()){
		while(all_results_.target_next()){
			pspline_.add_data(all_results_.get().get_param(),all_results_.get().get_S()->get_energy().get_x());
		}
		pso_info_.title("New PSpline run",'-');
		pspline_.compute_weights(); 

		double f;
		double f_tmp;
		bool is_min;
		Vector<double> param;
		Vector<unsigned int> idx;
		std::vector<Vector<double> > list_min;
#pragma omp parallel for private(is_min,idx,param,f,f_tmp)
		for(unsigned int i=0;i<parameter_space_.size();i++){
			is_min = true;
			idx = parameter_idx_[i];
			param = parameter_space_[i];
			f = pspline_.extrapolate(param);
			if(!isnan(f)){
				for(unsigned int j(0);j<Nfreedom_;j++){
					if(is_min){
						if(idx(j)+1<x_[j].size()){
							param(j) = x_[j](idx(j)+1);
							f_tmp = pspline_.extrapolate(param);
							if(isnan(f_tmp) || f>f_tmp){ is_min = false; }
						}
						if(is_min && idx(j)>0){
							param(j) = x_[j](idx(j)-1);
							f_tmp = pspline_.extrapolate(param);
							if(isnan(f_tmp) || f>f_tmp){ is_min = false; }
						}
						param(j) = x_[j](parameter_idx_[i](j));
					}
				}
			} else { is_min = false; }

			if(is_min){ 
#pragma omp critical(list_min)
				{
					list_min.push_back(param);
				}
			}
		}
		std::cout<<"min"<<std::endl;
		for(unsigned int i=0;i<list_min.size();i++){
			std::cout<<list_min[i]<<std::endl;
		}
#pragma omp parallel for
		for(unsigned int i=0;i<list_min.size();i++){
			bool tmp_test;
			std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(list_min[i]));
#pragma omp critical(all_results_)
			{
				if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
					sim->copy_S(all_results_.get().get_S()); 
					tmp_test=true;
				} else {
					sim->create_S(&system_param_);
					tmp_test=false;
				}
			}
			if(sim->is_created()){
				sim->run(tmp_test?10:1e6,tmax_);
#pragma omp critical(all_results_)
				{
					if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){
						all_results_.fuse_with_target(sim,MCSim::fuse);
					} else {
						all_results_.add_after_target(sim);
					}
				}
			}
		}
	}
}

bool VMCSpline::set_parameter_space(Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, bool const& border){
	Vector<double> tmp(Nfreedom_);
	for(unsigned int j(0); j<Nfreedom_;j++){ tmp(j) = x_[j](idx(j)); }
#pragma omp critical(parameter_space)
	{
		parameter_space_.push_back(tmp);
		parameter_idx_.push_back(idx);
	}
	if(border){
		for(unsigned int i(0);i<Nfreedom_;i++){
			if(i<Nfreedom_ && (idx(i) == 0 || idx(i) == x_[i].size()-1) ){
				i = Nfreedom_;
#pragma omp critical(border)
				{
					border_.push_back(tmp);
				}
			}
		}
	}

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<Nfreedom_;i++){
			if(idx(i) == x_[i].size()){ 
				if(i+1 == Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	} else {
		if(idx(0) == max0){ 
			if(1 == Nfreedom_){ return false; }
			idx(0) = min0;
			idx(1)++;
		}
		for(unsigned int i(1);i<Nfreedom_;i++){
			if(idx(i) == x_[i].size()){ 
				if(i+1 == Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return set_parameter_space(idx,min0,max0,border);
}

void VMCSpline::plot(){
	IOFiles data(get_filename()+".dat",true);
	all_results_.set_target();
	while(all_results_.target_next()){
		data<<all_results_.get().get_param()<<" "<<all_results_.get().get_S()->get_energy().get_x()<<IOFiles::endl;
	}

	IOFiles out(get_filename()+"-spline.dat",true);
	for(unsigned int i(0);i<parameter_space_.size();i++){
		out<<parameter_space_[i]<<" "<<pspline_.extrapolate(parameter_space_[i])<<IOFiles::endl;
	}

	Gnuplot plot("./",get_filename());
	plot.range("x","-2","2");
	plot.range("y","-2","2");
	plot.range("z","-1","");
	plot+="set ticslevel 0";
	plot+="splot '"+get_filename()+".dat' u 1:2:3 notitle,\\";
	plot+="      '"+get_filename()+"-spline.dat' u 1:2:3 t 'spline'";
	plot.save_file();
}
