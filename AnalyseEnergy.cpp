#include "AnalyseEnergy.hpp"

AnalyseEnergy::AnalyseEnergy(std::string const& sim):
	Analyse(sim)
{
	do_analyse();
}

void AnalyseEnergy::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); }
	if(level_==3){ 
		data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); 
		data_write_->precision(10);
	}
}

void AnalyseEnergy::close_files(){
	if(jd_write_){ 
		if(level_==7){ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); }
		rst_file_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

std::string AnalyseEnergy::extract_level_6(){ return extract_default(); }

std::string AnalyseEnergy::extract_level_5(){ return extract_default(); }

std::string AnalyseEnergy::extract_level_4(){ return extract_default(); }

std::string AnalyseEnergy::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	//(*read_)>>nof_;
	//jd_write_->add_to_header("\n");
	//jd_write_->write("number of jdfiles",nof_);
	//jd_write_->add_to_header("\n");
//

	CreateSystem cs(read_);
	cs.init(read_,this);
	cs.save();
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseEnergy::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$n^{-1}$'";
	gp.xrange(0,"");
	gp+="set key bottom";
	gp+="f(x)=a*x+b";
	gp+="a=1.0";
	gp+="b=1.0";
	gp+="set fit quiet";
	gp+="fit [0:0.01] f(x)  '"+filename_+".dat' u (1/$1):2:3 zerror via a,b";
	gp+="plot '"+filename_+".dat' u (1/$1):2:3 w e notitle,\\";
	gp+="     f(x) t sprintf('$E=%f$',b)";
	gp.save_file();
	gp.create_image(true);
	return filename_;
}

std::string AnalyseEnergy::extract_default(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	//jd_write_->add_to_header("\n");
	//jd_write_->write("number of jdfiles",nof_);
	//jd_write_->add_to_header("\n");
//
	CreateSystem cs(read_);
	cs.init(read_,this);
	cs.save();
	Data<double> E_;
	(*read_)>>E_;
	jd_write_->write("energy per site",E_);

	delete read_;
	read_ = NULL;

	return filename_;
}
