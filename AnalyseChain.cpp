#include "AnalyseChain.hpp"

AnalyseChain::AnalyseChain(std::string const& sim):
	Analyse(sim)
{}

void AnalyseChain::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); }
	if(level_==3 || level_==7){ data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); }
}

void AnalyseChain::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 7:{ 
					   if(nof_>1){/*if there is only one E data, there is no need to make a plot*/
						   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","E.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); 
					   }
				   } break;
			case 3:{ 
					   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+"-energy.png","energy.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+"-energy.gp",1000); 
					   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+"-polymerization.png","polymerization.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+"-polymerization.gp",1000); 
				   } break;
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

std::string AnalyseChain::extract_level_6(){
	/*E(param)|n=fixé*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	cs.init(read_,this);
	if(!all_link_names_.size()){ jd_write_->write("number of jdfiles",nof_); }
	jd_write_->add_to_header("\n");
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseChain::extract_level_5(){
	/*compare wavefunctions*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	unsigned int idx(0);
	double polymerization_strength;
	Data<double> E;
	Data<double> min_E;
	E.set_x(10.0);
	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E>>polymerization_strength;
		if(E.get_x()<min_E.get_x()){ 
			idx = i;
			min_E = E;
		}
	}
	delete read_;
	read_ = NULL;

	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E>>polymerization_strength;

		if(i==idx){
			cs.save(); 
			jd_write_->write("energy per site",min_E);
			jd_write_->write("polymerization strength",polymerization_strength);
		}
	}
	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseChain::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	double polymerization_strength;
	Data<double> E;

	CreateSystem cs(read_);
	cs.init(read_,this);
	(*read_)>>E>>polymerization_strength;

	cs.save(); 
	jd_write_->write("energy per site",E);
	jd_write_->write("polymerization strength",polymerization_strength);

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseChain::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	cs.init(read_,this);
	std::string link_name(cs.analyse(level_));
	cs.save();
	jd_write_->add_to_header("\n");

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseChain::extract_level_2(){
	Gnuplot gppolym(analyse_+path_+dir_,filename_+"-polymerization");
	gppolym+="set xlabel '$n^{-1}$'";
	gppolym.xrange(0,"");
	gppolym.yrange(0,"");
	gppolym+="set key bottom";
	gppolym+="plot '"+filename_+".dat' u (1/$1):6 t 'd-merization strength'";
	gppolym.save_file();
	gppolym.create_image(true);

	Gnuplot gpenergy(analyse_+path_+dir_,filename_+"-energy");
	gpenergy+="set xlabel '$n^{-1}$'";
	gpenergy.xrange(0,"");
	gpenergy+="set key bottom";
	gpenergy+="f(x)=a*x+b";
	gpenergy+="a=1.0";
	gpenergy+="b=1.0";
	gpenergy+="set fit quiet";
	gpenergy+="fit f(x) '"+filename_+".dat' u (1/$1):2:3 via a,b";
	gpenergy+="plot '"+filename_+".dat' u (1/$1):2:3 w e notitle,\\";
	gpenergy+="     f(x) t sprintf('$E=%f$',b)";
	gpenergy.save_file();
	gpenergy.create_image(true);
	return filename_;
}
