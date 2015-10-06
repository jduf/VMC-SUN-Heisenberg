#include "AnalyseMagnetization.hpp"

AnalyseMagnetization::AnalyseMagnetization(std::string const& sim, unsigned int const& max_level):
	Analyse(sim,max_level)
{}

void AnalyseMagnetization::open_files(){
	if(level_>1){ 
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
	}
	if(level_ == 4){
		jd_write_->write("number of jdfiles",nof_);
		jd_write_->add_header()->nl();
	}
	if(level_ == 6 || level_==4 || level_== 3){
		data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);
	}
}

void AnalyseMagnetization::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 4:{ list_rst_.last().figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","change_le_nom.png",RST::target(analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); } break;
			case 3:{ list_rst_.last().figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","change_le_nom.png",RST::target(analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000")); } break;
		}
		list_rst_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

std::string AnalyseMagnetization::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	System s(*read_);
	CreateSystem cs(&s);
	Vector<double> tmp(read_->read<Vector<double> >());
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);
	if(!all_link_names_.size()){ jd_write_->write("number of jdfiles",nof_); }
	jd_write_->add_header()->nl();
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

std::string AnalyseMagnetization::extract_level_5(){
	/*compare wavefunctions*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	unsigned int idx(0);
	Data<double> E;
	Data<double> min_E;
	min_E.set_x(10.0);
	Vector<unsigned int> M;
	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		Vector<double> tmp(read_->read<Vector<double> >());
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);
		(*read_)>>E;
		if(E.get_x()<min_E.get_x()){ 
			idx = i;
			min_E = E;
		}
	}
	delete read_;
	read_ = NULL;

	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	jd_write_->write("number of jdfiles",nof_);
	jd_write_->add_header()->nl();
	for(unsigned int i(0);i<nof_;i++){
		System s(*read_);
		CreateSystem cs(&s);
		Vector<double> tmp(read_->read<Vector<double> >());
		cs.init(&tmp,NULL);
		cs.set_IOSystem(this);
		(*read_)>>E;
		if(i==idx){
			s.save_input(*jd_write_); 
			jd_write_->write("energy per site",E);
		}
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMagnetization::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	Vector<double> tmp(read_->read<Vector<double> >());
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);
	s.save_input(*jd_write_);
	std::string link_name(cs.analyse(level_));
	jd_write_->add_header()->nl();

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMagnetization::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$M$' offset 0,1";
	gp+="set y2label '$\\dfrac{E}{n}$' rotate by 0";
	gp+="plot '"+filename_+".dat' u 1:($4==510?$2:1/0):3 w e t 'Fermi',\\";
	gp+="     '"+filename_+".dat' u 1:($4==511?$2:1/0):3 w e t 'Dirac',\\";
	gp+="     '"+filename_+".dat' u 1:($4==520?$2:1/0):3 w e t 'VBC'";
	gp.save_file();
	gp.create_image(true,true);

	Vector<unsigned int> ref;
	Vector<unsigned int> M;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	Data<double> E;
	Vector<unsigned int> vM(nof_);
	Vector<double> vE(nof_);
	vM(0) = 0;
	vE(0) = 2;
	for(unsigned int i(1);i<nof_;i++){
		(*read_)>>ref>>N>>m>>n>>M>>bc>>E;
		vE(i) = E.get_x();
		vM(i) = M(0);
	}

	Vector<unsigned int> index;
	vM.sort(std::less_equal<unsigned int>(),index);
	vE = vE.order(index);
	for(unsigned int i(1);i<nof_;i++){
		(*data_write_)<<-(vE(i)-vE(i-1))/(vM(i)-vM(i-1))<<" "<<1-vM(i)/double(vM(nof_-1))<<IOFiles::endl;
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string AnalyseMagnetization::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set key left";
	gp+="plot '"+filename_+".dat' u 1:2 notitle,\\";
	gp+="     1.0/3.0, 5.0/9.0, 7.0/9.0";
	gp.save_file();
	gp.create_image(true,true);
	return filename_;
}
