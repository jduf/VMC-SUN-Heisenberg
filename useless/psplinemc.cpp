/*! @file psplinemc.cpp */

#include "VMCSpline.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	VMCMinimization m(P);
	VMCSpline s(m);
	s.run(1);
	//List<MCSim> all_results_;
	//unsigned int Nfreedom(2);
	//double dx(0.04);
	//double min(-2.0);
	//double max(2.0);
	//Vector<double> param;
	//unsigned int tmax(P.get<unsigned int>("tmax"));
	//unsigned int N((max-min)/dx);
	//param.set(Nfreedom,min);
	//for(unsigned int dir(0);dir<Nfreedom;dir++){
//#pragma omp parallel for firstprivate(param)
		//for(unsigned int i=0;i<N;i+=N/10){
			//param(dir) = min+i*dx;
			//std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
			//sim->create_S(&P);
			//if(sim->is_created()){
				//sim->run(1e6,tmax);
//#pragma omp critical(all_results_)
				//{
					//if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){
						//all_results_.fuse_with_target(sim,MCSim::fuse);
					//} else {
						//all_results_.add_after_target(sim);
					//}
				//}
			//}
		//}
		//param(dir) = min;
	//}
	//param.set(Nfreedom,max);
	//for(unsigned int dir(0);dir<Nfreedom;dir++){
//#pragma omp parallel for firstprivate(param)
		//for(unsigned int i=0;i<N;i+=N/10){
			//param(dir) = min+i*dx;
			//std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
			//sim->create_S(&P);
			//if(sim->is_created()){
				//sim->run(1e6,tmax);
//#pragma omp critical(all_results_)
				//{
					//if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){
						//all_results_.fuse_with_target(sim,MCSim::fuse);
					//} else {
						//all_results_.add_after_target(sim);
					//}
				//}
			//}
		//}
		//param(dir) = max;
	//}
	//param.set();
	//for(unsigned int iter(0);iter<30;iter++){
//#pragma omp parallel for firstprivate(param)
		//for(unsigned int i=0;i<800;i++){
			//if(!param.ptr()){
				//param.set(Nfreedom);
				//set_param(param,min,Nfreedom,dx);
			//}
			//bool tmp_test;
			//std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
//
//#pragma omp critical(all_results_)
			//{
				//if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
					//sim->copy_S(all_results_.get().get_S()); 
					//tmp_test=true;
				//} else {
					//sim->create_S(&P);
					//tmp_test=false;
				//}
			//}
			//if(sim->is_created()){
				//sim->run(tmp_test?10:1e6,tmax);
//
//#pragma omp critical(all_results_)
				//{
					//if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){
						//all_results_.fuse_with_target(sim,MCSim::fuse);
					//} else {
						//all_results_.add_after_target(sim);
					//}
				//}
			//}
//
			//Matrix<double> c;
			//Vector<double> y;
//#pragma omp critical(all_results_)
			//{
				//c.set(all_results_.size(),Nfreedom);
				//y.set(c.row());
				//Vector<double> param_tmp;
				//unsigned int j(0);
				//while(all_results_.target_next()){
					//param_tmp = all_results_.get().get_param();
					//for(unsigned int k(0);k<Nfreedom;k++){
						//c(j,k) = param_tmp(k);
					//}
					//y(j) = all_results_.get().get_S()->get_energy().get_x();
					//j++;
				//}
			//}
			//if(sim->is_created() && !tmp_test){
				//PSpline s(2,c,y);
				//s.compute_weights();
				//unsigned int dir_opt(2*Nfreedom);
				//double old(sim->get_S()->get_energy().get_x());
				//double tmp;
				//for(unsigned int j(0);j<Nfreedom;j++){
					//param(j) += dx;
					//tmp = s.extrapolate(param);
					//if(s.extrapolate(param)<old){ 
						//dir_opt = 2*j; 
						//old = tmp;
					//}
					//param(j) -= 2*dx;
					//tmp = s.extrapolate(param);
					//if(tmp<old){ 
						//dir_opt = 2*j+1; 
						//old = tmp;
					//}
					//param(j) += dx;
				//}
				//if(dir_opt != 2*Nfreedom){
					//if(dir_opt%2){
						//dir_opt = (dir_opt-1)/2; 
						//param(dir_opt) -= dx;
					//} else {
						//dir_opt /= 2; 
						//param(dir_opt) += dx;
					//}
//
					//if(param(dir_opt) < min || param(dir_opt) > max){ set_param(param,min,Nfreedom,dx); } 
				//} else {
					//Rand<unsigned int> rnd(0,Nfreedom-1);
					//param(rnd.get()) += (2*rnd.get()>=Nfreedom?dx:-dx);
				//}
			//} else {
				//set_param(param,min,Nfreedom,dx);
			//}
		//}
		//param.set();
	//}
//
	//Matrix<double> c;
	//Vector<double> y;
	//c.set(all_results_.size(),Nfreedom);
	//y.set(c.row());
	//Vector<double> param_tmp;
	//unsigned int j(0);
	//IOFiles data("data.dat",true);
	//all_results_.set_target();
	//while(all_results_.target_next()){
		//param_tmp = all_results_.get().get_param();
		//for(unsigned int k(0);k<Nfreedom;k++){
			//c(j,k) = param_tmp(k);
		//}
		//y(j) = all_results_.get().get_S()->get_energy().get_x();
		//data<<param_tmp<<y(j)<<IOFiles::endl;
		//j++;
	//}
	//PSpline s(4,c,y);
	//s.compute_weights();
//
	//IOFiles out("attempt.dat",true);
	//Vector<double> tmp(Nfreedom);
	//for(unsigned int i(0);i<100;i++){
		//tmp(0) = min+i*dx;
		//for(unsigned int j(0);j<100;j++){
			//tmp(1) = min+j*dx;
			//out<<tmp<<" "<<s.extrapolate(tmp)<<IOFiles::endl;
		//}
	//}
//
	//Gnuplot plot("./","plot");
	//plot.range("x",-2,2);
	//plot.range("y",-2,2);
	//plot.range("z",-1,0);
	//plot+="set ticslevel 0";
	//plot+="splot 'data.dat' u 1:2:3 notitle,\\";
	//plot+="      'attempt.dat' u 1:2:3 notitle";
	//plot.save_file();
}

//void set_param(Vector<double>& param, double min, unsigned int Nfreedom, double dx){
	//Rand<unsigned int> rnd(1,100);
	//for(unsigned int j(0);j<Nfreedom;j++){
		//param(j) = min + rnd.get()*dx;
	//}
	//std::cerr<<"set param"<<std::endl;
//}
