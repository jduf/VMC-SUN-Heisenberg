#include "AnalyseSystematicLadder.hpp"

AnalyseSystematicLadder::AnalyseSystematicLadder(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd, bool const& display_results):
	Analyse(sim,path,max_level,run_cmd),
	display_results_(display_results),
	date_(Time().date("-"))
{
	do_analyse();
}

void AnalyseSystematicLadder::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);
	}
	switch(level_){
		case 6:
			{
				data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
				data_write_->precision(10);
			}break;
	}
}

void AnalyseSystematicLadder::close_files(){
	switch(level_){
		case 9:
			{
				best_wf_.set_target();
				best_sim_.set_target();
				Vector<double> param;
				double E(0);
				std::string sim("");
				while(best_wf_.target_next() && best_sim_.target_next() ){
					//std::cout<<best_sim_.get()<<" "<<best_wf_.get().get_param()<<" "<<best_wf_.get().get_energy()<<std::endl;
					if(E>best_wf_.get().get_energy().get_x()){
						E = best_wf_.get().get_energy().get_x();
						param = best_wf_.get().get_param();
						sim = best_sim_.get();
					}
				}
				best_wf_.set();
				best_sim_.set();
			}break;
		case 7:
			{
				best_J_.add_sort(kept_samples_J_.first_ptr(),MCSim::sort_by_theta_for_ladder);

				kept_samples_J_.set_target();
				unsigned int i(0);
				while(kept_samples_J_.target_next() && display_results_){
					i++;
					kept_samples_J_.get().display_results(my::tostring(i)+"-",sim_,info_,analyse_,path_,dir_,&list_rst_.last());
				}
				kept_samples_J_.set();
			}break;
		case 6:
			{
				kept_samples_all_.set_target();
				while(kept_samples_all_.target_next()){ kept_samples_all_.get().analyse(level_,this); }
				kept_samples_all_.set();

				std::string fname(dir_.substr(0,dir_.size()-1));
				Gnuplot gp(analyse_+path_,fname);
				gp.label("x","$\\theta$");
				gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
				gp.key("outside horizontal");
				gp+="plot '"+fname+".dat' u ($12==2?$5:1/0):6 pt 7 t 'Dimer $0$',\\";
				gp+="     '"+fname+".dat' u ($12==3?$5:1/0):6 pt 7 t 'Dimer $\\pi$',\\";
				gp+="     '"+fname+".dat' u ($12==4?$5:1/0):6 pt 7 t 'Square $00$',\\";
				gp+="     '"+fname+".dat' u ($12==5?$5:1/0):6 pt 7 t 'Square $0\\pi$',\\";
				gp+="     '"+fname+".dat' u ($12==6?$5:1/0):6 pt 7 t 'Square $\\pi\\pi$',\\";
				gp+="     '"+fname+".dat' u ($12==7?$5:1/0):6 pt 7 t 'Rectangle $\\pi\\pi\\pi$',\\";
				gp+="     '"+fname+".dat' u ($12==8?$5:1/0):6 pt 7 t 'Rectangle $0\\pi\\pi$',\\";
				gp+="     '"+fname+".dat' u ($12==9?$5:1/0):6 pt 7 t 'Rectangle $00\\pi$',\\";
				gp+="     '"+fname+".dat' u ($12==10?$5:1/0):6 pt 7 t 'Rectangle $000$'";
				gp.save_file();
				gp.create_image(true,"png");
				list_rst_.last().figure(rel_level_+analyse_+path_+fname+".png","Energy per site",RST::target(rel_level_+analyse_+path_+fname+".gp")+RST::width("1000"));

				best_J_.set_target();
				while(best_J_.target_next() && display_results_){ best_J_.get().display_results("",sim_,info_,analyse_,path_,dir_,&list_rst_.last()); }
				best_J_.set();
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
std::string AnalyseSystematicLadder::extract_level_9(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(read_->get_header());
	rst.save(true,false,true);

	List<MCSim> local_minima;
	VMCSystematic min(*read_);
	min.analyse(analyse_+path_+dir_,filename_,local_minima);

	local_minima.set_target();
	while(local_minima.target_next()){
		local_minima.get().save(*jd_write_);
		kept_samples_all_.add_end(local_minima.get_ptr());
		kept_samples_J_.add_sort(local_minima.get_ptr(),MCSim::sort_by_E);
		best_wf_.add_end(local_minima.get_ptr());
		best_sim_.add_end(std::make_shared<std::string>(read_->get_filename()));
	}

	delete read_;
	read_ = NULL;

	return filename_;
}
