#include "AnalyseExtract.hpp"

AnalyseExtract::AnalyseExtract(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file, unsigned int const& display_results, unsigned int const& ref):
	Analyse(sim,max_level,bash_file,ref),
	display_results_(display_results)
{
	if(display_results_){ std::cout<<"will call display_results() for "<<display_results_<<" samples"<<std::endl; }
	do_analyse();
}

AnalyseExtract::~AnalyseExtract(){
	if(study_){
		Gnuplot gp(analyse_+path_+dir_,sim_.substr(0,sim_.size()-1));
		gp.label("x","$\\frac{ 1}{N}$");
		gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
		gp+="plot '"+sim_.substr(0,sim_.size()-1)+".dat' u ($1/$2==6?1.0/$1:1/0):3 t '$k=6$',\\";
		gp+="     '"+sim_.substr(0,sim_.size()-1)+".dat' u ($1/$2==3?1.0/$1:1/0):3 t '$k=3$'";
		gp.save_file();
		gp.create_image(true,"png");
	}
}

void AnalyseExtract::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);
	}
	switch(level_){
		case 9:
			{
				jd_write_->write("number of different phase spaces",nof_);
				jd_write_->add_to_header()->np();
			}break;
		case 3:
			{
				jd_write_->write("number of different size",nof_);

				data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
				data_write_->precision(10);
			}break;
	}
}

void AnalyseExtract::close_files(){
	switch(level_){
		case 4:
			{ kept_samples_.set(); }break;
		case 3:
			{ list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","energy evolution with the system size",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); }break;
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
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(read_->get_header());
	rst.save(true,false,true);

	VMCExtract min(*read_,1e3,1e4);
	List<MCSim>::Node* target(min.analyse(analyse_+path_+dir_,filename_,kept_samples_));
	if(!kept_samples_.size()){ std::cerr<<__PRETTY_FUNCTION__<<" : need to have at least one sample in 'kept_samples_'"<<std::endl; }
	list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png",filename_,RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000"));
	if(target){ target->get()->save(*jd_write_); }

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseExtract::extract_level_8(){ return extract_best_of_previous_level(); }

std::string AnalyseExtract::extract_level_7(){ return extract_default(); }

std::string AnalyseExtract::extract_level_6(){
	kept_samples_.set_target();
	unsigned int i(0);
	while(kept_samples_.target_next() && i<display_results_){
		if(!i++){ std::cout<<std::string(6+path_.size()+dir_.size()+filename_.size(),' ')<<"|-> create lattice"; }
		else  { std::cout<<" "<<kept_samples_.size()-i<<std::flush; }

		kept_samples_.get().display_results(my::tostring(i)+"-",sim_,info_,analyse_,path_,dir_,&list_rst_.last(),false);
	}
	std::cout<<std::endl;

	return filename_;
}

std::string AnalyseExtract::extract_level_4(){
	double E(0.0);
	List<MCSim>::Node* target(NULL);

	kept_samples_.set_target();
	while(kept_samples_.target_next()){
		if(E>kept_samples_.get().get_energy().get_x()){
			E = kept_samples_.get().get_energy().get_x();
			target = kept_samples_.get_target();
		}
	}
	if(target){ target->get()->save(*jd_write_); }

	return filename_;
}

std::string AnalyseExtract::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Vector<double> tmp(*read_);
	System s(*read_);
	s.save(*data_write_);

	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);

	jd_write_->add_to_header()->nl();
	cs.save(*jd_write_);

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseExtract::extract_level_2(){ return fit_therm_limit(); }
