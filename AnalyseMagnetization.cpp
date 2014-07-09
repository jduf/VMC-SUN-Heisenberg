#include "AnalyseMagnetization.hpp"

void AnalyseMagnetization::open_files(){
	if(level_>1){ 
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
	}
	if(level_ == 6 || level_==4){
		data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);
	}
}

void AnalyseMagnetization::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 4:{ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","change_le_nom.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
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

std::string AnalyseMagnetization::extract_level_7(){
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

std::string AnalyseMagnetization::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	cs.init(read_,this);
	if(!all_link_names_.size()){ (*jd_write_)("number of jdfiles",nof_); }
	jd_write_->add_to_header("\n");
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseMagnetization::extract_level_5(){
	/*compare wavefunctions*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	unsigned int idx(0);
	Data<double> E;
	Data<double> min_E;
	min_E.set_x(10.0);
	Vector<unsigned int> M;
	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E;
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
	jd_write_->add_to_header("\n");
	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E;
		if(i==idx){
			cs.save(); 
			(*jd_write_)("energy per site",E);
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMagnetization::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	CreateSystem cs(read_);
	cs.init(read_,this);
	cs.save();
	std::string link_name(cs.analyse(level_));
	jd_write_->add_to_header("\n");

	return filename_;
}

std::string AnalyseMagnetization::extract_level_3(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$M$' offset 0,1";
	gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
	gp+="plot '"+filename_+".dat' u 1:($4==510?$2:1/0):3 w e t 'Fermi',\\";
	gp+="     '"+filename_+".dat' u 1:($4==511?$2:1/0):3 w e t 'Dirac',\\";
	gp+="     '"+filename_+".dat' u 1:($4==520?$2:1/0):3 w e t 'VBC'";
	gp.save_file();
	gp.create_image(true);

	return filename_;
}
