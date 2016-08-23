#include "AnalyseHoneycomb.hpp"

AnalyseHoneycomb::AnalyseHoneycomb(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd):
	Analyse(sim,path,max_level,run_cmd)
{ do_analyse(); }

void AnalyseHoneycomb::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false); }
	if(level_ == 7){
		data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
		data_write_->precision(10);
	}
}

void AnalyseHoneycomb::close_files(){
	if(jd_write_){
		if(level_==7){ list_rst_.last().figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Honeycomb",RST::target(analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); }
		list_rst_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

std::string AnalyseHoneycomb::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Vector<double> tmp(*read_);
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}
