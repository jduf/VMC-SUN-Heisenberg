#include "AnalyseLadder.hpp"

AnalyseLadder::AnalyseLadder(std::string const& path, unsigned int const& max_level):
	Analyse(path,max_level)
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
		case 9:
			{
				rst_file_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png","Parameter sets",RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000")); 
			} break;
		case 6:
			{
				rst_file_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); 
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
	IOFiles in(sim_+path_+dir_+filename_+".jdbin",false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(in.get_header()); 
	rst.save(false,true); 

	VMCMinimization min(in);
	min.find_save_and_plot_minima(10,*jd_write_,analyse_+path_+dir_,filename_);

	return filename_;
}

std::string AnalyseLadder::extract_level_8(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	Vector<double> tmp;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>tmp;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);

		/*!Draw the lattice with the witdth related to t_*/
		cs.lattice(info_+path_+dir_,filename_+"-"+my::tostring(i));
		rst_file_.last().figure(dir_+filename_+"-"+my::tostring(i)+".png",RST::math("E="+my::tostring(s.get_energy().get_x())+"\\pm"+my::tostring(s.get_energy().get_dx())),RST::target(dir_+filename_+"-"+my::tostring(i)+".pdf")+RST::scale("200")); 

		if(!i){/*only the best set of parameter is kept*/
			cs.save_param(*jd_write_);
			s.save_input(*jd_write_);
			s.save_output(*jd_write_);
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
	Vector<double> tmp;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>tmp;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);

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
		(*read_)>>tmp;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);
		if(i==idx){
			cs.save_param(*jd_write_);
			s.save_input(*jd_write_);
			s.save_output(*jd_write_);
			i=nof_;
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	Vector<double> tmp(read_->read<Vector<double> >());
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);
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
	gp.key("left Left");
	gp+="plot '"+filename_+".dat' u 5:6";
	gp.save_file();
	gp.create_image(true,true);

	return filename_;
}
