#include "AnalyseSystematic.hpp"

AnalyseSystematic::AnalyseSystematic(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file, unsigned int const& ref):
	Analyse(sim,max_level,bash_file,ref)
{
	do_analyse();
}

void AnalyseSystematic::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);
	}
	switch(level_){
		case 9:
			{
				jd_write_->write("number of different wavefunction",nof_);
				jd_write_->add_to_header()->np();
			}break;
		case 3:
			{
				data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
				data_write_->precision(10);

				jd_write_->write("all minima",kept_samples_.size());
				kept_samples_.set_target();
				while(kept_samples_.target_next()){
					kept_samples_.get().save(*jd_write_);
					kept_samples_.get().save(*data_write_);
				}
			}break;
	}
}

void AnalyseSystematic::close_files(){
	switch(level_){
		case 3:
			{ list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","energy evolution with the system size",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); }break;
		case 2:
			{
				kept_samples_.first().analyse(level_,this);
				kept_samples_.set();
			}break;
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
std::string AnalyseSystematic::extract_level_9(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(read_->get_header());
	rst.save(true,false,true);

	List<MCSim> local_minima;
	VMCSystematic vmc(*read_);
	vmc.analyse(local_minima);
	vmc.plot(analyse_+path_+dir_,filename_);
	list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png",filename_,RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000"));

	local_minima.set_target();
	while(local_minima.target_next()){
		local_minima.get().save(*jd_write_);
		kept_samples_.add_sort(local_minima.get_ptr(),MCSim::sort_by_E);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}
