#include "AnalyseExtract.hpp"

AnalyseExtract::AnalyseExtract(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd):
	Analyse(sim,path,max_level,run_cmd)
{
	do_analyse();
}

void AnalyseExtract::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
	}

	//if(level_==8){
	//jd_write_->write("number of different number of dof",nof_);
	//jd_write_->add_header()->np();
	//}

	if(level_==3){
		data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);
		data_write_->precision(10);
	}
}

void AnalyseExtract::close_files(){
	switch(level_){
		case 4:
			{ kept_samples_.set(); }
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
std::string AnalyseExtract::extract_level_9(){
	IOFiles in(sim_+path_+dir_+filename_+".jdbin",false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(in.get_header());
	rst.save(false,true);

	VMCExtract min(in);
	min.plot(analyse_+path_+dir_,filename_,kept_samples_);
	list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png",filename_,RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000")); 

	return filename_;
}

std::string AnalyseExtract::extract_level_8(){
	IOFiles in(sim_+path_+dir_+filename_+".jdbin",false);

	kept_samples_.set_target();
	while(kept_samples_.target_next()){
		kept_samples_.get().display_results(sim_,info_,analyse_,path_,dir_,&list_rst_.last());
	}
	//rst.text(iof.get_header());

	return filename_;
}

std::string AnalyseExtract::extract_level_4(){
	kept_samples_.set_target();
	double E(0.0);
	List<MCSim>::Node* target(NULL);
	while(kept_samples_.target_next()){
		if(E>kept_samples_.get().get_energy().get_x()){
			E = kept_samples_.get().get_energy().get_x();
			target = kept_samples_.get_target();
		}
	}
	if(target){ target->get()->write(*jd_write_); }

	return filename_;
}

std::string AnalyseExtract::extract_level_3(){
	IOFiles in(sim_+path_+dir_+filename_+".jdbin",false);
	MCSim sim(in);
	sim.save(*data_write_);


	return filename_;
}

std::string AnalyseExtract::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.range("x","0","");
	gp.label("x","$\\dfrac{1}{n}$");
	gp.label("y2","$\\dfrac{E}{nN^2}$","rotate by 0");
	gp.key("left Left");
	gp+="plot '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) w e notitle";
	gp.save_file();
	gp.create_image(true,true);

	list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png","energy",RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000")); 
	return filename_;
}
