#include "AnalyseMin.hpp"

AnalyseMin::AnalyseMin(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd):
	Analyse(sim,path,max_level,run_cmd)
{
	child_in_AnalyseMin_ = true;
	do_analyse();
}

void AnalyseMin::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);
	}
	switch(level_){
		case 8:
			{
				jd_write_->write("number of different wavefunction",nof_);
				jd_write_->add_header()->np();
			}break;
		case 7:
			{
				jd_write_->write("number of different wavefunction",nof_);
				jd_write_->add_header()->np();
			}break;
		case 5:
			{ jd_write_->write("number of different boundary condition",nof_); }break;
		case 3:
			{
				jd_write_->write("number of different size",nof_);

				data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
				data_write_->precision(10);
			}break;
	}
}

void AnalyseMin::close_files(){
	switch(level_){
		case 9:
			{ list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png","Parameter sets",RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000")); } break;
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

std::string AnalyseMin::extract_level_9(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(read_->get_header());
	rst.save(false,true);

	VMCMinimization min(*read_,true,"ANA");
	min.find_save_and_plot_minima(10,*jd_write_,analyse_+path_+dir_,filename_);

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMin::extract_level_8(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	(*read_)>>nof_;

	Vector<double> tmp;
	std::string tmp_filename(filename_);
	for(unsigned int i(0);i<nof_;i++){
		if(!i){ std::cout<<std::string(6+path_.size()+dir_.size()+filename_.size(),' ')<<"|-> create lattice"; }
		else  { std::cout<<" "<<nof_-i<<std::flush; }

		(*read_)>>tmp;
		System s(*read_);
		s.create_cluster(true);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		rst_file_ = &list_rst_.last();
		filename_ = tmp_filename+"-"+my::tostring(i);
		cs.set_IOSystem(this);
		cs.display_results();
		rst_file_ = NULL;

		/*only the best (first) set of parameter is kept*/
		if(!i){ cs.save(*jd_write_); }
	}
	std::cout<<std::endl;

	filename_ = tmp_filename;

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMin::extract_level_7(){ return extract_best_of_previous_level(); }

std::string AnalyseMin::extract_level_6(){ return extract_best_of_previous_level(); }

std::string AnalyseMin::extract_level_5(){ return extract_default(); }

std::string AnalyseMin::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	(*read_)>>nof_;
	/*!must save now nof_ because it doesn't refer to the number of file in
	 * the next directory but in the next-next directory*/
	jd_write_->write("number of different boundary condition",nof_);

	for(unsigned int i(0);i<nof_;i++){
		Vector<double> tmp(*read_);
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);

		jd_write_->add_header()->nl();
		cs.save(*jd_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMin::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		Vector<double> tmp(*read_);
		System s(*read_);
		s.save(*data_write_);

		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);

		jd_write_->add_header()->nl();
		cs.save(*jd_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMin::extract_level_2(){ return fit_therm_limit(); }
