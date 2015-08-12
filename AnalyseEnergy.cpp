#include "AnalyseEnergy.hpp"

AnalyseEnergy::AnalyseEnergy(std::string const& path):
	Analyse(path)
{
	do_analyse();
}

void AnalyseEnergy::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); 
		if(level_==5){ jd_write_->write("number of different boundary condition",nof_); }
		if(level_==3 || level_==7){
			data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); 
			data_write_->precision(10);
		}
	}
}

void AnalyseEnergy::close_files(){
	if(jd_write_){ 
		if(level_==7){ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); }
		if(level_==3){ rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); }
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

/*different wavefunction*/
std::string AnalyseEnergy::extract_level_5(){ return extract_default(); }

/*different boundary condition*/
std::string AnalyseEnergy::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	/*!must save now nof_ because it doesn't refer to the number of file in the
	 * next directory but in the next-next directory*/
	jd_write_->write("number of different boundary condition",nof_); 

	Data<double> E;
	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		cs.construct_GenericSystem(read_,this);
		(*read_)>>E;

		jd_write_->add_header()->nl();
		s.save_input(*jd_write_); 
		jd_write_->write("energy per site",E);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseEnergy::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		cs.construct_GenericSystem(read_,this);
		std::string link_name(cs.analyse(level_));

		jd_write_->add_header()->nl();
		s.save_input(*jd_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseEnergy::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$n^{-1}$'";
	gp.range("x","0","");
	gp+="set key bottom";
	gp+="f(x)=a*x+b";
	gp+="a=1.0";
	gp+="b=1.0";
	gp+="set fit quiet";
	gp+="fit [0:0.01] f(x)  '"+filename_+".dat' u (1/$1):2:3 zerror via a,b";
	gp+="plot '"+filename_+".dat' u (1/$1):($6==1?$2:1/0):3 w e t 'P',\\";
	gp+="     '"+filename_+".dat' u (1/$1):($6==-1?$2:1/0):3 w e t 'A',\\";
	gp+="     '"+filename_+".dat' u (1/$1):($6==0?$2:1/0):3 w e t '0',\\";
	gp+="     f(x) t sprintf('$E=%f$',b)";
	gp.save_file();
	gp.create_image(true,true);
	return filename_;
}

std::string AnalyseEnergy::extract_default(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	jd_write_->add_header()->nl();
	Data<double> E;
	System s(*read_);
	CreateSystem cs(&s);
	cs.construct_GenericSystem(read_,this);
	if(level_ == 6){ (*read_)>>nof_; }
	(*read_)>>E;

	s.save_input(*jd_write_);
	jd_write_->write("energy per site",E);

	delete read_;
	read_ = NULL;

	return filename_;
}
