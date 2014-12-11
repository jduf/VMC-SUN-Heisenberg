#include "AnalyseHoneycomb.hpp"

AnalyseHoneycomb::AnalyseHoneycomb(std::string const& sim):
	Analyse(sim)
{
	do_analyse();
}

void AnalyseHoneycomb::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); }
	if(level_ == 7){ 
		data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); 
		data_write_->precision(10);
	}
}

void AnalyseHoneycomb::close_files(){
	if(jd_write_){ 
		if(level_==7){ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Honeycomb",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); }
		rst_file_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

std::string AnalyseHoneycomb::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	cs.init(read_,this);
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}
