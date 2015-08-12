#include "AnalyseLadder.hpp"

AnalyseLadder::AnalyseLadder(std::string const& path):
	Analyse(path)
{
	do_analyse();
}

void AnalyseLadder::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); 

		if(level_==8){ 
			jd_write_->write("yahhooo",nof_); 
			jd_write_->add_header()->np();
		} 

		if(level_==6){ 
			data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); 
			data_write_->precision(10);
		} 
	}
}

void AnalyseLadder::close_files(){
	switch(level_){
		case 6:
			{
				rst_file_.last().link_figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","allks",rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); 
			} break;
	}
	if(jd_write_){ 
		rst_file_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

std::string AnalyseLadder::extract_level_9(){
	VMCMinimization min(sim_+path_+dir_+filename_+".jdbin");
	nof_ = 5;
	jd_write_->write("number of compared simulations",nof_);
	min.save_best(nof_,*jd_write_);
	/*will also plot the E(param)*/

	return filename_;
}

std::string AnalyseLadder::extract_level_8(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		cs.construct_GenericSystem(read_,this);

		if(!i){
			s.save_input(*jd_write_);
			s.save_output(*jd_write_);
			cs.save_param(*jd_write_);
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_7(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	double E_tmp(0);
	unsigned int idx(0);
	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		cs.construct_GenericSystem(read_,this);

		if(E_tmp>s.get_energy().get_x()){
			E_tmp = s.get_energy().get_x();
			idx=i;
		}
	}

	delete read_;
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	Data<double> E;
	(*read_)>>nof_;
	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		cs.construct_GenericSystem(read_,this);
		if(i==idx){
			s.save_input(*jd_write_);
			s.save_output(*jd_write_);
			cs.save_param(*jd_write_);
			i=nof_;
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	System s(*read_);
	CreateSystem cs(&s);
	cs.construct_GenericSystem(read_,this);
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_5(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.range("x","0","pi/2");
	gp.label("x","$\\theta$");
	gp.label("y2","$\\dfrac{E}{n}$","rotate by 0");
	gp+="plot '"+filename_+".dat' u 5:6";
	gp.save_file();
	gp.create_image(true,true);

	return filename_;
}
