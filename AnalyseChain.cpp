#include "AnalyseChain.hpp"

AnalyseChain::AnalyseChain(std::string const& sim):
	Analyse(sim)
{
	std::cout<<"Will proceed to the analyse of a chain system. It"<<std::endl;
	std::cout<<"will consist of an analyse of :"<<std::endl;
	std::cout<<"+ first neighbour correlations"<<std::endl;
	std::cout<<"+ long range correlations"<<std::endl;
	std::cout<<"+ extraction of the energy"<<std::endl;
	std::cout<<"The energies will be compared and the best"<<std::endl;
	std::cout<<"simulation will be kept. There will the be a"<<std::endl;
	std::cout<<"comparison of the energies and the polymerization"<<std::endl;
	std::cout<<"strength in function of the system size."<<std::endl;
}

void AnalyseChain::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); 
		if(level_==6){ jd_write_->write("number of wavefunction",nof_); } 
		if(level_==5){ jd_write_->write("number of boundary condition",nof_); }
		if(level_==3 || level_==7){ 
			data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true); 
			data_write_->precision(10);
		}
		if(level_==3){ (*data_write_)<<"%N m bc n E(x,dx,#,conv) d-strength exponents"<<IOFiles::endl; }
	}
}

void AnalyseChain::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 7:{ 
					   if(nof_>1){/*if there is only one E data, there is no need to make a plot*/
						   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy per site",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); 
					   }
				   } break;
			case 3:{ 
					   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+"-energy.png","Energy per site",analyse_+path_+dir_.substr(0,dir_.size()-1)+"-energy.gp",1000); 
					   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+"-exponents.png","Critical Exponents",analyse_+path_+dir_.substr(0,dir_.size()-1)+"-exponents.gp",1000); 
					   rst_file_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+"-polymerization.png","Polymerization",analyse_+path_+dir_.substr(0,dir_.size()-1)+"-polymerization.gp",1000); 
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

/*calls cs.analyse(level) to plot E(delta)*/
std::string AnalyseChain::extract_level_6(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	CreateSystem cs(read_);
	cs.init(read_,this);
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

/*different wavefunction*/
std::string AnalyseChain::extract_level_5(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	unsigned int idx(0);
	double polymerization_strength;
	Vector<double> exponents;
	Data<double> min_E;
	Data<double> E;
	min_E.set_x(10.0);
	double diff_e(0);
	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E>>polymerization_strength>>exponents;
		diff_e = std::abs(E.get_x() - min_E.get_x());
		if(E.get_x()<min_E.get_x()){ 
			idx = i;
			min_E = E;
		}
	}
	if(nof_==2){
		jd_write_->add_to_header("\n");
		jd_write_->add_to_header("difference of energy between the wavefunction = " +tostring(diff_e));
	}
	delete read_;
	read_ = NULL;

	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E>>polymerization_strength>>exponents;

		if(i==idx){
			jd_write_->add_to_header("\n");
			cs.save(); 
			jd_write_->write("energy per site",min_E);
			jd_write_->write("polymerization strength",polymerization_strength);
			jd_write_->write("critical exponents",exponents);
		}
	}
	delete read_;
	read_ = NULL;

	return filename_;
}

/*different boundary condition*/
std::string AnalyseChain::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	jd_write_->write("number of boundary condition",nof_); 

	double polymerization_strength;
	Vector<double> exponents;
	Data<double> E;

	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E>>polymerization_strength>>exponents;

		jd_write_->add_to_header("\n");
		cs.save(); 
		jd_write_->write("energy per site",E);
		jd_write_->write("polymerization strength",polymerization_strength);
		jd_write_->write("critical exponents",exponents);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*different magnetization*/
std::string AnalyseChain::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;

	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		std::string link_name(cs.analyse(level_));

		jd_write_->add_to_header("\n");
		cs.save();
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*different system size*/
std::string AnalyseChain::extract_level_2(){
	Gnuplot gpenergy(analyse_+path_+dir_,filename_+"-energy");
	gpenergy+="set xlabel '$n^{-1}$'";
	gpenergy+="set y2label '$\\dfrac{E}{n}$' rotate by 0";
	gpenergy+="set key bottom";
	gpenergy.xrange(0,"");
	gpenergy+="f(x)=a*x**b+E0";
	gpenergy+="a=-1.0";
	gpenergy+="b=2.0";
	gpenergy+="E0=-1.0";
	gpenergy+="set fit quiet";
	gpenergy+="fit f(x) '"+filename_+".dat' u (1/$4):($3==1?$5:1/0):6 zerror via a,b,E0";
	gpenergy+="plot '"+filename_+".dat' u (1/$4):($3==1?$5:1/0):6 w e t 'P',\\";
	gpenergy+="     '"+filename_+".dat' u (1/$4):($3==0?$5:1/0):6 w e t 'O',\\";
	gpenergy+="     f(x) t sprintf('$E=%f$',E0)";
	gpenergy.save_file();
	gpenergy.create_image(true);

	Gnuplot gppolym(analyse_+path_+dir_,filename_+"-polymerization");
	gppolym+="set xlabel '$n^{-1}$'";
	gppolym+="set y2label 'd-merization strength'";
	gppolym+="set key bottom";
	gppolym.xrange(0,"");
	gppolym.yrange(0,"");
	gppolym+="plot '"+filename_+".dat' u (1/$4):($3==1?$9:1/0) t 'P',\\";
	gppolym+="     '"+filename_+".dat' u (1/$4):($3==0?$9:1/0) t 'O'";
	gppolym.save_file();
	gppolym.create_image(true);

	Gnuplot gpexp(analyse_+path_+dir_,filename_+"-exponents");
	gpexp+="set lmargin at screen 0.15";
	gpexp+="set rmargin at screen 0.85";
	gpexp+="";
	gpexp+="set multiplot";
	gpexp+="set tmargin at screen 0.55";
	gpexp+="set bmargin at screen 0.15";
	gpexp+="set xrange [0:] writeback";
	gpexp+="set xlabel '$n^{-1}$'";
	gpexp+="set xtics nomirror";
	gpexp+="set y2tics";
	gpexp+="unset ytics";
	gpexp+="plot '"+filename_+".dat' u (1/$4):($3==1?$13:1/0) lc 1 t '$a_0=2$'";
	gpexp+="";
	gpexp+="set tmargin at screen 0.95";
	gpexp+="set bmargin at screen 0.55";
	gpexp+="set xrange restore";
	gpexp+="set ytics";
	gpexp+="unset xtics";
	gpexp+="unset y2tics";
	gpexp+="unset xlabel";
	gpexp+="plot '"+filename_+".dat' u (1/$4):($3==1?$11:1/0) lc 2 t '$a_{\\pi}=2-2/N$'";
	gpexp+="unset multiplot";
	gpexp.save_file();
	gpexp.create_image(true);

	return filename_;
}
