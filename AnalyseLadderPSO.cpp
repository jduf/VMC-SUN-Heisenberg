#include "AnalyseLadderPSO.hpp"

AnalyseLadderPSO::AnalyseLadderPSO(std::string const& sim, unsigned int const& max_level, bool const& run_cmd, unsigned int const& display_results):
	Analyse(sim,max_level,run_cmd,2),
	display_results_(display_results)
{ do_analyse(); }

void AnalyseLadderPSO::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);

		if(level_==8){
			jd_write_->write("number of different dof",nof_);
			jd_write_->add_to_header()->np();
		}

		if(level_==6){
			data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
			data_write_->precision(10);
		}
	}
}

void AnalyseLadderPSO::close_files(){
	switch(level_){
		case 9:
			{ list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png","Parameter sets",RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000")); } break;
		case 6:
			{ list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); } break;
	}
	if(jd_write_){
		list_rst_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

/*extract VMCMinimization and plot*/
std::string AnalyseLadderPSO::extract_level_9(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	List<MCSim> local_minima;
	VMCExtract(*read_,0,0).analyse(analyse_+path_+dir_,filename_,local_minima);

	local_minima.set_target();
	unsigned int i(0);
	/*Need to save at least one file*/
	(*jd_write_)<<std::min(local_minima.size(),display_results_?display_results_:1);
	while(local_minima.target_next() && i++<display_results_+1){
		local_minima.get().save(*jd_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*show best solutions*/
std::string AnalyseLadderPSO::extract_level_8(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	rst_file_ = &list_rst_.last();
	(*read_)>>nof_;

	Vector<double> tmp;
	std::string tmp_filename(filename_);
	for(unsigned int i(0);i<nof_;i++){
		/*!Draw the lattice with the witdth related to t_*/
		if(!i){ std::cout<<std::string(6+path_.size()+dir_.size()+filename_.size(),' ')<<"|-> create lattice"; }
		else { std::cout<<" "<<nof_-i<<std::flush; }

		(*read_)>>tmp;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		filename_ = tmp_filename+"-"+my::tostring(i);
		cs.set_IOSystem(this);
		cs.display_results();

		if(!i){/*only the best set of parameter is kept*/
			cs.save(*jd_write_);
		}
	}
	std::cout<<std::endl;

	filename_ = tmp_filename;

	delete read_;
	read_ = NULL;
	rst_file_ = NULL;

	return filename_;
}

/*compare wavefunction (different ref_)*/
std::string AnalyseLadderPSO::extract_level_7(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
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
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Data<double> E;
	(*read_)>>nof_;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>tmp;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);
		if(i==idx){
			cs.save(*jd_write_);
			i=nof_;
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*plot energy for different J*/
std::string AnalyseLadderPSO::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Vector<double> tmp(read_->read<Vector<double> >());
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	rst_file_ = &list_rst_.last();
	replace_title_with_link_in_rst_ = true;
	cs.set_IOSystem(this);
	std::string link_name(cs.analyse(level_));
	rst_file_ = NULL;
	replace_title_with_link_in_rst_ = false;

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseLadderPSO::extract_level_5(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.range("x","0","pi/2");
	gp.label("x","$\\theta$");
	gp.label("y2","$\\dfrac{E}{n}$","rotate by 0");
	gp.key("left Left");
	gp+="plot '"+filename_+".dat' u 5:6";
	gp.save_file();
	gp.create_image(true,"png");

	return filename_;
}
