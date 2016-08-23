#include "AnalyseChain.hpp"

AnalyseChain::AnalyseChain(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd):
	Analyse(sim,path,max_level,run_cmd),
	outfile_(NULL)
{
	std::cout<<"Will proceed to the analyse SU(N) chains. It will consist "
		"of an analyse of :"<<std::endl<<
		"+ first neighbour correlations"<<std::endl<<
		"+ long range correlations"<<std::endl<<
		"+ extraction of the energy"<<std::endl<<
		"The energies will be compared and the best simulation will be "
		"kept. There will the be a"<<std::endl<<
		"comparison of the energies and the polymerization strength in "
		"function of the system size."<<std::endl;

	if(path==""){ outfile_ = new IOFiles("sun-chains.dat",true,false); }
	do_analyse();
}

AnalyseChain::~AnalyseChain(){
	if(outfile_){ delete outfile_; }
}

void AnalyseChain::open_files(){
	if(level_>1){
		jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true,false);
		if(level_==7){
			jd_write_->write("number of different wavefunction",nof_);
			jd_write_->add_header()->np();
		}
		if(level_==6){
			jd_write_->write("number of different boundary condition",nof_);
			jd_write_->add_header()->np();
		}
		if(level_==3 || level_==8){
			data_write_ = new IOFiles(analyse_+path_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
			data_write_->precision(10);
		}
		if(level_==3){ (*data_write_)<<"%N m bc n E(x,dx,#,conv) d-strength exponents"<<IOFiles::endl; }
	}
}

void AnalyseChain::close_files(){
	if(jd_write_){
		switch(level_){
			case 8:
				{
					if(nof_>1){/*if there is only one E data, there is no need to make a plot*/
						list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","Energy per site",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp")+RST::width("1000"));
					}
				} break;
			case 3:
				{
					list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+"-energy.png","Energy per site",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+"-energy.gp")+RST::width("1000"));
					list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+"-exponents.png","Critical Exponents",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+"-exponents.gp")+RST::width("1000"));
					list_rst_.last().figure(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+"-polymerization.png","Polymerization",RST::target(rel_level_+analyse_+path_+dir_.substr(0,dir_.size()-1)+"-polymerization.gp")+RST::width("1000"));
				} break;
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

/*{*//*!Extract a single simulation

	   Extracts the results and plots the correlations, long range correlations
	   and structure factor
	   *//*}*/
std::string AnalyseChain::extract_level_8(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Vector<double> tmp(*read_);
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);

	if(!all_link_names_.size()){
		jd_write_->add_header()->title("Simulations",'<');
		jd_write_->write("number of jdfiles",nof_);
		(*data_write_)<<"% t E dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	}
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}

/*calls cs.analyse(level) to plot E(ti)*/
std::string AnalyseChain::extract_level_7(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	unsigned int idx(0);
	std::string link_name(filename_);

	double E(666);
	Vector<double> t;

	double d_tmp;
	Vector<double> vd_tmp;

	(*read_)>>nof_;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>t;
		System s(*read_);
		(*read_)>>d_tmp>>vd_tmp;

		if(s.get_energy().get_x()<E){
			E = s.get_energy().get_x();
			idx = i;
		}
	}

	delete read_;
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	(*read_)>>nof_;
	for(unsigned int i(0);i<idx+1;i++){
		(*read_)>>t;
		System s(*read_);

		if(i==idx){
			CreateSystem cs(&s);
			cs.init(&t,NULL);
			cs.set_IOSystem(this);
			link_name = cs.analyse(level_);
		} else {
			(*read_)>>d_tmp>>vd_tmp;
		}
	}

	delete read_;
	read_ = NULL;

	return link_name;
}

/*different wavefunction*/
std::string AnalyseChain::extract_level_6(){
	std::string link_name;
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	Vector<unsigned int> ref;
	Vector<unsigned int> M;
	(*read_)>>nof_;

	switch(nof_){
		case 1:
			{
				double polymerization_strength;
				Vector<double> exponents;
				Vector<double> t;
				(*read_)>>t;

				System s(*read_);
				(*read_)>>polymerization_strength>>exponents;

				CreateSystem cs(&s);
				cs.init(&t,NULL);
				cs.set_IOSystem(this);

				cs.save(*jd_write_);
				jd_write_->write("polymerization strength",polymerization_strength);
				jd_write_->write("critical exponents",exponents);

				if(outfile_){
					s.save(*outfile_);
					(*outfile_)<<t<<" "<<exponents<<" "<<polymerization_strength<<IOFiles::endl;
					//(*outfile_)<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<ti<<" ";
					//(*outfile_)<<E<<" "<<exponents<<" "<<polymerization_strength<<IOFiles::endl;
				}

				std::cerr<<__PRETTY_FUNCTION__<<" only one wf"<<std::endl;
			}break;
		case 2:
			{
				/*!This part should load the two different wavefunctions. If
				 * the ChainPolymerization has a vector set to 1, there this
				 * wavefunction is equivalent to the ChainFermi one. In that
				 * case the energy obtained with this wavefunction is correct
				 * and should replace the one obtained by ChainFermi.
				 * ChainFermi should therefore be the wavefunction to save but
				 * with the energy of the ChainPolymerized one. In the other
				 * case, the energy given by ChainPolymerized should also be
				 * taken even if the one obtained by ChainFermi is lower (the
				 * Fermi one is then wrong) and ChainPolymerized is the
				 * wavefunction to save.
				 *
				 * As we know that SU(9) m=3 is gapless, it will save the
				 * correct critical exponent
				 unsigned int eta_idx(( (N==9&&m==3) || ref(2)==0)?0:1);
				 */
				std::cerr<<__PRETTY_FUNCTION__<<" : don't know what to do"<<std::endl;
			}break;
		default:{ std::cerr<<__PRETTY_FUNCTION__<<" : too many wavefunctions"<<std::endl; }
	}
	return filename_;
}

/*different J_*/
std::string AnalyseChain::extract_level_5(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	(*read_)>>nof_;
	Vector<double> t(*read_);
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&t,NULL);
	cs.set_IOSystem(this);

	cs.save(*jd_write_);
	jd_write_->write("polymerization strength",read_->read<double>());
	jd_write_->write("critical exponents",read_->read<Vector<double> >());

	return filename_;
}

/*different boundary condition*/
std::string AnalyseChain::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	/*!must save now nof_ because it doesn't refer to the number of file in the
	 * next directory but in the next-next directory*/
	jd_write_->write("number of different boundary condition",nof_);

	std::cout<<nof_<<std::endl;
	Vector<double> t;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>t;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&t,NULL);
		cs.set_IOSystem(this);

		cs.save(*jd_write_);
		jd_write_->write("polymerization strength",read_->read<double>());
		jd_write_->write("critical exponents",read_->read<Vector<double> >());
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*different magnetization*/
std::string AnalyseChain::extract_level_3(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
	(*read_)>>nof_;

	Vector<double> t;
	for(unsigned int i(0);i<nof_;i++){
		(*read_)>>t;
		System s(*read_);
		CreateSystem cs(&s);
		cs.init(&t,NULL);
		cs.set_IOSystem(this);

		cs.analyse(level_);

		s.save(*jd_write_);
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*different system size*/
std::string AnalyseChain::extract_level_2(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Vector<unsigned int> ref(*read_);
	unsigned int N(read_->read<unsigned int>());
	unsigned int m(read_->read<unsigned int>());

	delete read_;
	read_ = NULL;

	Gnuplot gpenergy(analyse_+path_+dir_,filename_+"-energy");
	gpenergy.label("x","$n^{-1}$");
	gpenergy.label("y2","$\\dfrac{E}{n}$","rotate by 0");
	gpenergy.range("x","0","");
	gpenergy+="set key bottom";
	if((N/m)%2 == 0){
		if(N/m!=2){
			gpenergy+="f(x)=a1*x**b1+E1";
			gpenergy+="a1=1.0";
			gpenergy+="b1=2.0";
			gpenergy+="E1=1.0";
			gpenergy+="g(x)=a2*x**b2+E2";
			gpenergy+="a2=-1.0";
			gpenergy+="b2=2.0";
			gpenergy+="E2=1.0";
			gpenergy+="set fit quiet";
			gpenergy+="fit f(x) '"+filename_+".dat' u (1/$4):(($3==1 && floor($4*$2/$1)%2==0)?$5:1/0):6 zerror via a1,b1,E1";
			gpenergy+="fit g(x) '"+filename_+".dat' u (1/$4):(($3==1 && floor($4*$2/$1)%2==1)?$5:1/0):6 zerror via a2,b2,E2";
			gpenergy+="plot '"+filename_+".dat' u (1/$4):($3==1?$5:1/0):6 w e t 'P',\\";
			gpenergy+="     '"+filename_+".dat' u (1/$4):($3==0?$5:1/0):6 w e t 'O',\\";
			gpenergy+="     f(x) lc 3 t sprintf('$E=%f$',E1),\\";
			gpenergy+="     g(x) lc 3 t sprintf('$E=%f$',E2)";
		} else {
			gpenergy+="plot '"+filename_+".dat' u (1/$4):($3==1?$5:1/0):6 w e t 'P',\\";
			gpenergy+="     '"+filename_+".dat' u (1/$4):($3==0?$5:1/0):6 w e t 'O'";
		}
	} else {
		gpenergy+="f(x)=a*x**b+E";
		gpenergy+="a=1.0";
		gpenergy+="b=2.0";
		gpenergy+="E=1.0";
		gpenergy+="set fit quiet";
		gpenergy+="fit f(x) '"+filename_+".dat' u (1/$4):($3==1?$5:1/0):6 zerror via a,b,E";
		gpenergy+="plot '"+filename_+".dat' u (1/$4):($3==1?$5:1/0):6 w e t 'P',\\";
		gpenergy+="     '"+filename_+".dat' u (1/$4):($3==0?$5:1/0):6 w e t 'O',\\";
		gpenergy+="     f(x) lc 3 t sprintf('$E=%f$',E)";
	}
	gpenergy.save_file();
	gpenergy.create_image(true,true);

	Gnuplot gpexp(analyse_+path_+dir_,filename_+"-exponents");
	gpexp.multiplot();
	gpexp.margin("0.15","0.85","0.95","0.55");
	gpexp+="unset xtics";
	gpexp.range("x","0","0.037");
	gpexp.range("y","1.3","1.82");
	gpexp+="plot '"+filename_+".dat' u (1/$4):($3==1?$11:1/0) lc 1 t '$\\nu$',\\";
	gpexp+="     "+my::tostring(2.0-2.0/N)+" lc 1 notitle";
	gpexp+="";
	gpexp.margin("0.15","0.85","0.55","0.15");
	gpexp+="set xtics nomirror";
	gpexp.label("x","$n^{-1}$");
	gpexp+="unset ytics";
	gpexp.range("y");
	gpexp.range("y2","1.5","2.5");
	gpexp+="set y2tics mirror";
	gpexp+="plot '"+filename_+".dat' u (1/$4):($3==1?$13:1/0) axes x1y2 lc 2 t '$\\nu$',\\";
	gpexp+="     2.0 axes x1y2 lc 2 notitle";
	gpexp+="unset multiplot";
	gpexp.save_file();
	gpexp.create_image(true,true);

	Gnuplot gppolym(analyse_+path_+dir_,filename_+"-polymerization");
	gppolym.label("x","$n^{-1}$");
	gppolym.label("y2","d-merization strength");
	gppolym.range("x",0,"");
	gppolym.range("y",0,"");
	gppolym+="set key bottom";
	gppolym+="plot '"+filename_+".dat' u (1/$4):($3==1?$9:1/0) t 'P',\\";
	gppolym+="     '"+filename_+".dat' u (1/$4):($3==0?$9:1/0) t 'O'";
	gppolym.save_file();
	gppolym.create_image(true,true);

	return filename_;
}
