#include "VMCSpline.hpp"

VMCSpline::VMCSpline(Parseur& P):
	VMCMinimization(P),
	pspline_(2)
{}

void VMCSpline::compute_border(){
	unsigned int N;
	Vector<double> param(Nfreedom_);
	for(unsigned int dir(0);dir<Nfreedom_;dir++){
		N = x_[dir].size();
		std::cout<<dir<<" "<<x_[dir]<<std::endl;
/*{Description*/
/*!
		all permutation of 
			[0,0,z]
			[0,1,z]
			[1,0,z]
			[1,1,z]
			[0,y,0]
			[0,y,1]
			[1,y,0]
			[1,y,1]
			[x,0,0]
			[x,0,1]
			[x,1,0]
			[x,1,1]
 */
/*}*/

#pragma omp parallel for firstprivate(param)
		for(unsigned int i=0;i<N;i+=N/10){
			param(dir) = x_[dir](i);
#pragma omp critical
			{
			std::cout<<param<<std::endl;
			}
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
}

void VMCSpline::run(){
	set_time();
	pso_info_.title("New PSpline run",'-');

	if(all_results_.size()){
		compute_border();
		pspline_.compute_weights(); 
		Vector<double> param;
		Vector<unsigned int> idx;
		std::vector<Vector<double> > list_min;
#pragma omp parallel for
			for(unsigned int i=0;i<x_[0].size();i++){
				split vector x_ and check all points and find all min
			}
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

#pragma omp critical(add_to_vector_)
			{
				pspline_.add_data(list_min[i],sim->get_S()->get_energy().get_x());
			}

			if(sim->is_created() && !tmp_test){
				if(!omp_get_thread_num()){
					std::cout<<"recompute weights"<<std::endl;
					pspline_.compute_weights(); 
				}
				unsigned int dir_opt(2*Nfreedom_);
				double old(sim->get_S()->get_energy().get_x());
				double tmp;
				for(unsigned int dir(0);dir<Nfreedom_;dir++){
					if(idx(dir)+1<x_[dir].size()){
						param(dir) = x_[dir](idx(dir)+1);
						tmp = pspline_.extrapolate(param);
						if(pspline_.extrapolate(param)<old){ 
							dir_opt = 2*dir; 
							old = tmp;
						}
					}
					if(idx(dir)>0){
						param(dir) = x_[dir](idx(dir)-1);
						tmp = pspline_.extrapolate(param);
						if(tmp<old){ 
							dir_opt = 2*dir+1; 
							old = tmp;
						}
					}
					param(dir) = x_[dir](idx(dir));
				}

#pragma omp critical(add_results_)
				{
					std::cout<<"do"<<dir_opt<<std::endl;
				}
				if(dir_opt != 2*Nfreedom_){
					if(dir_opt%2){
						dir_opt = (dir_opt-1)/2; 
						idx(dir_opt)--;
					} else {
						dir_opt /= 2; 
						idx(dir_opt)++;
					}
					param(dir_opt) = x_[dir_opt](idx(dir_opt));

					if(param(dir_opt) < x_[dir_opt](0) || param(dir_opt) > x_[dir_opt](x_[dir_opt].size())){ set_param(param,idx); } 
				} else {
					Rand<unsigned int> rnd(0,Nfreedom_-1);
					dir_opt = rnd.get();
					idx(dir_opt) += (2*rnd.get()>=Nfreedom_?1:-1);
					param(dir_opt) = x_[dir_opt](idx(dir_opt));
				}
			} else {
				set_param(param,idx);
			}
		}
	}
}

void VMCSpline::set_param(Vector<double>& param, Vector<unsigned int>& idx){
	for(unsigned int dir(0);dir<Nfreedom_;dir++){
		Rand<unsigned int> rnd(0,x_[dir].size()-1);
		idx(dir) = rnd.get();
		param(dir) = x_[dir](idx(dir));
	}
	std::cerr<<"set param"<<std::endl;
}

void VMCSpline::plot(){
	IOFiles data("data.dat",true);
	all_results_.set_target();
	while(all_results_.target_next()){
		data<<all_results_.get().get_param()<<all_results_.get().get_S()->get_energy().get_x()<<IOFiles::endl;
	}

	pspline_.compute_weights();
	IOFiles out("attempt.dat",true);
	Vector<double> tmp(Nfreedom_);
	for(unsigned int i(0);i<x_[0].size();i++){
		tmp(0) = x_[0](i);
		for(unsigned int j(0);j<x_[1].size();j++){
			tmp(1) = x_[1](j);
			out<<tmp<<" "<<pspline_.extrapolate(tmp)<<IOFiles::endl;
		}
	}

	Gnuplot plot("./","plot");
	plot.range("x",-2,2);
	plot.range("y",-2,2);
	plot.range("z",-1,0);
	plot+="set ticslevel 0";
	plot+="splot 'data.dat' u 1:2:3 notitle,\\";
	plot+="      'attempt.dat' u 1:2:3 notitle";
	plot.save_file();
}
