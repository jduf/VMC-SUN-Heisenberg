#include "AnalyseEnergy.hpp"

AnalyseEnergy::AnalyseEnergy(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd):
	Analyse(sim,path,max_level,run_cmd)
{
	do_analyse();
}

AnalyseEnergy::~AnalyseEnergy(){
	Gnuplot gp(analyse_+path_+dir_,sim_.substr(0,sim_.size()-1));
	gp.label("x","$\\frac{1}{N}$");
	gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
	gp+="plot '"+sim_.substr(0,sim_.size()-1)+".dat' u ($1/$2==6?1.0/$1:1/0):3 t '$k=6$',\\";
	gp+="     '"+sim_.substr(0,sim_.size()-1)+".dat' u ($1/$2==3?1.0/$1:1/0):3 t '$k=3$'";
	gp.save_file();
	gp.create_image(true,true);
}

void AnalyseEnergy::open_files(){
	if(level_){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
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
				data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);
				data_write_->precision(10);
			}break;
	}
}

void AnalyseEnergy::close_files(){
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

std::string AnalyseEnergy::extract_level_8(){ return extract_default(); }

std::string AnalyseEnergy::extract_level_7(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int idx(0);
	double E(666);
	Vector<double> param;

	(*read_)>>nof_;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>param;
		System s(*read_);
		if(s.get_energy().get_x()<E){
			E = s.get_energy().get_x();
			idx = i;
		}
	}

	delete read_;
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	(*read_)>>nof_;
	for(unsigned int i(0);i<idx+1;i++){
		(*read_)>>param;
		System s(*read_);

		if(i==idx){
			CreateSystem cs(&s);
			cs.init(&param,NULL);
			cs.set_IOSystem(this);
			cs.save(*jd_write_);
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseEnergy::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int idx(0);
	double E(666);
	Vector<double> param;

	(*read_)>>nof_;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>param;
		System s(*read_);
		if(s.get_energy().get_x()<E){
			E = s.get_energy().get_x();
			idx = i;
		}
	}

	delete read_;
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	(*read_)>>nof_;
	for(unsigned int i(0);i<idx+1;i++){
		(*read_)>>param;
		System s(*read_);

		if(i==idx){
			CreateSystem cs(&s);
			cs.init(&param,NULL);
			cs.set_IOSystem(this);
			cs.save(*jd_write_);
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*different wavefunction*/
std::string AnalyseEnergy::extract_level_5(){ return extract_default(); }

/*different boundary condition*/
std::string AnalyseEnergy::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
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

std::string AnalyseEnergy::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		Vector<double> tmp(*read_);
		System s(*read_);
		s.save(*data_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseEnergy::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="f(x) = a*x+b";
	gp+="set fit quiet";
	gp+="fit [0:0.025] f(x) '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) yerror via a,b";
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

std::string AnalyseEnergy::extract_default(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	Vector<double> tmp(*read_);
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);
	if(level_ == 6){ (*read_)>>nof_; }

	jd_write_->add_header()->nl();
	cs.save(*jd_write_);

	delete read_;
	read_ = NULL;

	return filename_;
}
