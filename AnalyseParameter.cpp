#include "AnalyseParameter.hpp"


void AnalyseParameter::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); }
	if(level_==3 || level_==6){ data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); }
}

void AnalyseParameter::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 6:{ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","E.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
			case 3:{ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","d-merization.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
		}
		rst_file_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

std::string AnalyseParameter::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	/*might be a way to init and set_IOSystem at the same time, but first need
	 * to check if won't be a problem somewhere else*/
	cs.init(read_);
	cs.set_IOSystem(this);
	/*save only once the general datas*/
	/*will save twice delta...*/
	if(!all_link_names_.size()){ 
		cs.save(*jd_write_);
		jd_write_->add_to_header("\n");
		(*jd_write_)("number of jdfiles",nof_);
		jd_write_->add_to_header("\n");
	}
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseParameter::extract_level_5(){
	/*E(param)|n=fixé*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	/*might be a way to init and set_IOSystem at the same time, but first need
	 * to check if won't be a problem somewhere else*/
	cs.init(read_);
	cs.set_IOSystem(this);
	/*save only once the general datas*/
	/*will save twice delta...*/
	if(!all_link_names_.size()){ (*jd_write_)("number of jdfiles",nof_); }
	jd_write_->add_to_header("\n");
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseParameter::extract_level_4(){
	/*compare wavefunctions*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	unsigned int idx(0);
	double polymerization_strength;
	Data<double> E;
	Data<double> min_E;
	E.set_x(10.0);
	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		/*might be a way to init and set_IOSystem at the same time, but first need
		 * to check if won't be a problem somewhere else*/
		cs.init(read_);
		cs.set_IOSystem(this);
		(*read_)>>E>>polymerization_strength;
		if(E.get_x()<min_E.get_x()){ 
			idx = i;
			min_E = E;
		}
	}
	delete read_;
	read_ = NULL;

	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	(*jd_write_)("number of jdfiles",nof_);
	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		/*might be a way to init and set_IOSystem at the same time, but first need
		 * to check if won't be a problem somewhere else*/
		cs.init(read_);
		cs.set_IOSystem(this);
		(*read_)>>E>>polymerization_strength;

		if(i==idx){
			jd_write_->add_to_header("\n");
			cs.save(*jd_write_); 
			(*jd_write_)("energy per site",min_E);
			(*jd_write_)("polymerization strength",polymerization_strength);
		}
	}
	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseParameter::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	CreateSystem cs(read_);
	/*might be a way to init and set_IOSystem at the same time, but first need
	 * to check if won't be a problem somewhere else*/
	cs.init(read_);
	cs.set_IOSystem(this);
	(*jd_write_)("number of jdfiles",nof_);
	jd_write_->add_to_header("\n");
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseParameter::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$n^{-1}$'";
	gp.xrange(0,"");
	gp.yrange(0,"");
	gp+="set key bottom";
	gp+="plot '"+filename_+".dat' u (1/$1):2 t 'd-merization strength'";
	gp.save_file();
	gp.create_image(true);
	return filename_;
}

