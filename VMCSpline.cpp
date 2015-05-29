#include "VMCSpline.hpp"

VMCSpline::VMCSpline(Parseur& P):
	VMCMinimization(P,"PSpline"),
	pspline_(2),
	border_(NULL)
{}

VMCSpline::~VMCSpline(){
	if(border_){ delete[] border_; }
}

void VMCSpline::init(bool border){
	set_time();
	pso_info_.title("New PSpline run",'-');
	(void)(border);
	//if(border){
		//pso_info_.text("set border"+RST::nl_);
		//border_ = new Vector<double>[compute_border_size(n)];
//#pragma omp parallel
		//{
			//unsigned int min0(omp_get_thread_num()*x_[0].size()/omp_get_num_threads());
			//unsigned int max0((omp_get_thread_num()+1)*x_[0].size()/omp_get_num_threads());
			//Vector<unsigned int> idx;
			//idx.set(Nfreedom_,0);
			//idx(0) = min0;
			///*could remove the while I think*/
			//while(set_parameter_space(x_,idx,min0,max0,border));
		//}
		//unsigned int incr(border_.size()/12);
//#pragma omp parallel for
		//for(unsigned int i=0;i<border_.size();i+=incr){
			//compute_vmc(border_[i]);
		//}
	//}
}

void VMCSpline::run(){
	if(all_results_.size()){
		while(all_results_.target_next()){
			pspline_.add_data(all_results_.get().get_param(),all_results_.get().get_S()->get_energy().get_x());
		}
		pso_info_.title("New PSpline run",'-');
		pspline_.compute_weights(); 

		std::cout<<"try to find min"<<std::endl;
#pragma omp parallel 
		{
			unsigned int min0(omp_get_thread_num()*x_[0].size()/omp_get_num_threads());
			unsigned int max0((omp_get_thread_num()+1)*x_[0].size()/omp_get_num_threads());
			Vector<unsigned int> idx;
			idx.set(Nfreedom_,0);
			idx(0) = min0;
			go_through_parameter_space(x_,idx,min0,max0,&VMCSpline::run_if_min);
		}
		for(unsigned int j(0); j<Nfreedom_;j++){ 
			std::cout<<min_idx_[j]<<std::endl;
		}

#pragma omp parallel for
		for(unsigned int i=0;i<min_idx_.size();i++){
			//Vector<double>* x(new Vector<double>[Nfreedom_]);
			//for(unsigned int j(0);j<Nfreedom_;j++){
				//x[j].set(5,0);
				//for(unsigned int k(0);k<x[j].size();k++){
					//if(min_idx_[i](j) >= k-2 && min_idx_[i](j)+k-2 < .size()){
						//x[j](k) = x_(min_idx_[i](j)+k-2);
					//}
				//}
			//}
			//Vector<unsigned int> idx;
//
			//std::vector<Vector<double> > parameter_space;
			//std::vector<Vector<unsigned int> > parameter_idx;
			//set_parameter_space(x,parameter_space,parameter_idx,idx);
			//for(unsigned int j(0);j<parameter_space.size();j++){
				//compute_vmc(parameter_space[j]);
			//}
			//delete[] x;

			Vector<double> param(Nfreedom_);
			for(unsigned int j(0); j<Nfreedom_;j++){ param(j) = x_[j](min_idx_[i](j)); }
			compute_vmc(param);
		}
	}
}

bool VMCSpline::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSpline::*f)(Vector<double>*, Vector<unsigned int> const&)){
	(this->*f)(x,idx);

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<Nfreedom_;i++){
			if(idx(i) == x[i].size()){ 
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
			if(idx(i) == x[i].size()){ 
				if(i+1 == Nfreedom_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return go_through_parameter_space(x,idx,min0,max0,f);
}

void VMCSpline::run_if_min(Vector<double>* x, Vector<unsigned int> const& idx){
	double f;
	double f_tmp;
	bool is_min(true);
	Vector<double> param(Nfreedom_);
	for(unsigned int i(0); i<Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	f = pspline_.extrapolate(param);
	if(!isnan(f)){
		for(unsigned int j(0);j<Nfreedom_;j++){
			if(is_min){
				if(idx(j)+1<x_[j].size()){
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
			min_idx_.push_back(idx);
		}
	}
}

void VMCSpline::save_spline_data(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(Nfreedom_);
	for(unsigned int i(0); i<Nfreedom_;i++){ param(i) = x[i](idx(i)); }
	(*out_)<<param<<" "<<pspline_.extrapolate(param)<<IOFiles::endl;
}

void VMCSpline::plot(){
	out_ = new IOFiles(get_filename()+".dat",true);
	all_results_.set_target();
	while(all_results_.target_next()){
		(*out_)<<all_results_.get().get_param()<<" "<<all_results_.get().get_S()->get_energy().get_x()<<IOFiles::endl;
	}
	delete out_;

	out_ = new IOFiles(get_filename()+"-spline.dat",true);
	Vector<unsigned int> idx(Nfreedom_,0);
	go_through_parameter_space(x_,idx,0,0,&VMCSpline::save_spline_data);
	delete out_;
	out_ = NULL;

	Gnuplot plot("./",get_filename());
	plot.range("x","-2","2");
	plot.range("y","-2","2");
	plot.range("z","-1","");
	plot+="set ticslevel 0";
	plot+="splot '"+get_filename()+".dat' u 1:2:3 notitle,\\";
	plot+="      '"+get_filename()+"-spline.dat' u 1:2:3 t 'spline'";
	plot.save_file();
}
