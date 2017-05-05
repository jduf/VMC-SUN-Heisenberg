#include "AnalyseLadder.hpp"

AnalyseLadder::AnalyseLadder(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file):
	Analyse(sim,max_level,bash_file)
{
	do_analyse();
}

AnalyseLadder::~AnalyseLadder(){
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

void AnalyseLadder::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);
	}
	switch(level_){
		case 8:
			{
				jd_write_->write("number of different wavefunction",nof_);
				jd_write_->add_to_header()->np();
			}break;
		case 7:
			{
				jd_write_->write("number of different wavefunction",nof_);
				jd_write_->add_to_header()->np();
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

void AnalyseLadder::close_files(){
	switch(level_){
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

std::string AnalyseLadder::extract_level_8(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	rst_file_ = list_rst_.last_ptr().get();

	Vector<double> tmp(*read_);
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);

	jd_write_->add_to_header()->nl();
	cs.save(*jd_write_);
	cs.display_results();

	delete read_;
	read_ = NULL;
	rst_file_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_7(){ return extract_best_of_previous_level(); }

std::string AnalyseLadder::extract_level_6(){ return extract_best_of_previous_level(); }

std::string AnalyseLadder::extract_level_5(){ return extract_default(); }

std::string AnalyseLadder::extract_level_4(){
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

		jd_write_->add_to_header()->nl();
		cs.save(*jd_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		Vector<double> tmp(*read_);
		System s(*read_);
		//s.save(*data_write_);

		CreateSystem cs(&s);
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);

		jd_write_->add_to_header()->nl();
		cs.save(*jd_write_);
		cs.save(*data_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseLadder::extract_level_2(){ return fit_therm_limit(); }
