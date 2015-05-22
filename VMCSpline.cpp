#include "VMCSpline.hpp"

VMCSpline::VMCSpline(Parseur& P):
	VMCMinimization(P),
	pspline_(2)
{
}

void VMCSpline::compute_border(double const& border, unsigned int const& dir){
	Vector<double> param(Nfreedom_,border);
	unsigned int N((max_-min_)/dx_);

#pragma omp parallel for firstprivate(param)
	for(unsigned int i=0;i<N;i+=N/10){
		param(dir) = min_+i*dx_;
		std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
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

void VMCSpline::set_param(Vector<double>& param){
	Rand<unsigned int> rnd(1,100);
	for(unsigned int j(0);j<Nfreedom_;j++){
		param(j) = min_ + rnd.get()*dx_;
	}
	std::cerr<<"set param"<<std::endl;
}

void VMCSpline::run(){
	for(unsigned int i(0);i<Nfreedom_;i++){
		compute_border(min_,i);
		compute_border(max_,i);
	}

	Vector<double> param;
	for(unsigned int iter(0);iter<30;iter++){
#pragma omp parallel for firstprivate(param)
		for(unsigned int i=0;i<800;i++){
			if(!param.ptr()){
				param.set(Nfreedom_);
				set_param(param);
			}
			bool tmp_test;
			std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));

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

#pragma omp critical(add_results_)
			{
				pspline_.add_data(param,sim->get_S()->get_energy().get_x());
			}

			if(sim->is_created() && !tmp_test){
				pspline_.compute_weights();
				unsigned int dir_opt(2*Nfreedom_);
				double old(sim->get_S()->get_energy().get_x());
				double tmp;
				for(unsigned int j(0);j<Nfreedom_;j++){
					param(j) += dx_;
					tmp = pspline_.extrapolate(param);
					if(pspline_.extrapolate(param)<old){ 
						dir_opt = 2*j; 
						old = tmp;
					}
					param(j) -= 2*dx_;
					tmp = pspline_.extrapolate(param);
					if(tmp<old){ 
						dir_opt = 2*j+1; 
						old = tmp;
					}
					param(j) += dx_;
				}
				if(dir_opt != 2*Nfreedom_){
					if(dir_opt%2){
						dir_opt = (dir_opt-1)/2; 
						param(dir_opt) -= dx_;
					} else {
						dir_opt /= 2; 
						param(dir_opt) += dx_;
					}

					if(param(dir_opt) < min_ || param(dir_opt) > max_){ set_param(param); } 
				} else {
					Rand<unsigned int> rnd(0,Nfreedom_-1);
					param(rnd.get()) += (2*rnd.get()>=Nfreedom_?dx_:-dx_);
				}
			} else {
				set_param(param);
			}
		}
		param.set();
	}
}
