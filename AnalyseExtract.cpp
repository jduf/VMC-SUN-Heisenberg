#include "AnalyseExtract.hpp"

AnalyseExtract::AnalyseExtract(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd):
	Analyse(sim,path,max_level,run_cmd)
{
	do_analyse();
}

AnalyseExtract::~AnalyseExtract(){
	Gnuplot gp(analyse_+path_+dir_,sim_.substr(0,sim_.size()-1));
	gp.label("x","$\\frac{1}{N}$");
	gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
	gp+="plot '"+sim_.substr(0,sim_.size()-1)+".dat' u ($1/$2==6?1.0/$1:1/0):3 t '$k=6$',\\";
	gp+="     '"+sim_.substr(0,sim_.size()-1)+".dat' u ($1/$2==3?1.0/$1:1/0):3 t '$k=3$'";
	gp.save_file();
	gp.create_image(true,true);
}

void AnalyseExtract::open_files(){
	if(level_){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
	}
	switch(level_){
		case 3:
			{
				data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);
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
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	RSTFile rst(info_+path_+dir_,filename_);
	rst.text(read_->get_header());
	rst.save(false,true);

	VMCExtract min(*read_,true);
	min.select_minima_and_plot(analyse_+path_+dir_,filename_,kept_samples_);
	list_rst_.last().figure(rel_level_+analyse_+path_+dir_+filename_+".png",filename_,RST::target(rel_level_+analyse_+path_+dir_+filename_+".gp")+RST::width("1000"));

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseExtract::extract_level_7(){
	kept_samples_.set_target();
	while(kept_samples_.target_next()){
		kept_samples_.get().display_results(sim_,info_,analyse_,path_,dir_,&list_rst_.last());
	}

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
	if(target){ target->get()->save(*jd_write_); }

	return filename_;
}

std::string AnalyseExtract::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	Vector<double> tmp(*read_);
	System s(*read_);
	s.save(*data_write_);

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseExtract::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="f(x) = a*x+b";
	gp+="set fit quiet";
	gp+="fit f(x) '"+filename_+".dat' u ($3!=72?1.0/$3:1/0):($5/($1*$1)):($6/($1*$1)) yerror via a,b";
	gp+="set print \"../"+sim_.substr(0,sim_.size()-1)+".dat\" append";
	gp+="print \"`head -1 '"+filename_+".dat' | awk '{print $1 \" \" $2}'`\",\" \",b";
	gp.range("x","0","");
	gp.label("x","$\\frac{1}{n}$");
	gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
	gp+="plot '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) w e notitle,\\";
	gp+="     f(x) t sprintf('%f',b)";
	gp.save_file();
	gp.create_image(true,true);

	return filename_;
}
