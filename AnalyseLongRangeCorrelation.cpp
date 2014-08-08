#include "AnalyseLongRangeCorrelation.hpp"

AnalyseLongRangeCorrelation::AnalyseLongRangeCorrelation(std::string const& sim):
	Analyse(sim)
{}

void AnalyseLongRangeCorrelation::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); }
	if(level_==7){ data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); }
}

void AnalyseLongRangeCorrelation::close_files(){
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

std::string AnalyseLongRangeCorrelation::extract_level_7(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	cs.init(read_,this);
	/*Only one call of cs.save() is needed*/
	if(!all_link_names_.size()){ 
		cs.save();
		jd_write_->add_to_header("\n");
		(*jd_write_)("number of jdfiles",nof_);
		jd_write_->add_to_header("\n");
	}
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}
