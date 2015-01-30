#include "AnalyseChain.hpp"

AnalyseChain::AnalyseChain(std::string const& path):
	Analyse(path)
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

	if(path==""){ outfile_ = new IOFiles("sun-chains.dat",true); }
	do_analyse();
}

AnalyseChain::~AnalyseChain(){
	if(outfile_){ delete outfile_; }
}

void AnalyseChain::open_files(){
	if(level_>1){ jd_write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true); 
		if(level_==6){ jd_write_->write("number of different wavefunction",nof_); } 
		if(level_==5){ jd_write_->write("number of different boundary condition",nof_); }
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

/*calls cs.analyse(level) to plot E(ti)*/
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
	std::string link_name;
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	Vector<unsigned int> ref;
	Vector<unsigned int> M;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;

	jd_write_->add_header()->nl();
	switch(nof_){
		case 1:
			{ 
				double polymerization_strength;
				Vector<double> exponents;
				Data<double> E;
				Vector<double> ti;
				(*read_)>>ref>>N>>m>>n>>M>>bc;
				if(ref(2) == 1){ (*read_)>>ti; }
				else{ ti.set(N/m,1); }
				(*read_)>>E>>polymerization_strength>>exponents;
				if(ref(2)==1){
					ChainPolymerized chain(ref,N,m,n,M,bc,ti);
					chain.set_IOSystem(this);
					chain.save();
				} else {
					switch(ref(1)){
						case 1:
							{
								ChainFermi<double> chain(ref,N,m,n,M,bc);
								chain.set_IOSystem(this);
								chain.save();
							}break;
						case 2:
							{
								ChainFermi<std::complex<double> > chain(ref,N,m,n,M,bc);
								chain.set_IOSystem(this);
								chain.save();
							}break;
						default:{ std::cerr<<"std::string AnalyseChain::extract_level_5() : ref undefined"<<std::endl; }
					}
				}
				if(outfile_){
					(*outfile_)<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<ti<<" ";
					(*outfile_)<<E<<" "<<exponents<<" "<<polymerization_strength<<IOFiles::endl;
				}
				jd_write_->write("energy per site",E);
				jd_write_->write("polymerization strength",polymerization_strength);
				jd_write_->write("critical exponents",exponents);
				std::cerr<<"only one wf for "<<N<<" "<<m<<" "<<n<<" "<<bc<<std::endl; 
			}break;
		case 2:
			{ 
				double polymerization_strength[2];
				Vector<double> exponents[2];
				Data<double> E[2];
				Vector<double> ti;
				for(unsigned int i(0);i<nof_;i++){
					(*read_)>>ref>>N>>m>>n>>M>>bc;
					if(ref(2) == 1){ (*read_)>>ti; }
					else{ ti.set(N/m,1); }
					(*read_)>>E[ref(2)]>>polymerization_strength[ref(2)]>>exponents[ref(2)];
				}

				std::cerr<<"find condition when poly==fermi"<<std::endl;
				if(!are_equal(ti(N/m-1),1)){
					ref(1) = 1;
					ref(2) = 1;
					ChainPolymerized chain(ref,N,m,n,M,bc,ti);
					chain.set_IOSystem(this);
					chain.save();
				} else {
					ref(2) = 0;
					switch(ref(1)){
						case 1:
							{
								ChainFermi<double> chain(ref,N,m,n,M,bc);
								chain.set_IOSystem(this);
								chain.save();
							}break;
						case 2:
							{
								ChainFermi<std::complex<double> > chain(ref,N,m,n,M,bc);
								chain.set_IOSystem(this);
								chain.save();
							}break;
						default:{ std::cerr<<"std::string AnalyseChain::extract_level_5() : ref undefined"<<std::endl; }
					}
				}

				/*As we know that SU(9) m=3 is gapless, it will save the correct critical exponent*/
				std::cerr<<"find condition when poly==fermi"<<std::endl;
				unsigned int eta_idx(((N==9&&m==3) || are_equal(ti(N/m-1),1))?0:1);
				if(outfile_){
					(*outfile_)<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<ti<<" ";
					(*outfile_)<<E[1]<<" "<<exponents[eta_idx]<<" "<<polymerization_strength[ref(2)]<<IOFiles::endl;
				}
				jd_write_->write("energy per site",E[1]);
				jd_write_->write("polymerization strength",polymerization_strength[ref(2)]);
				jd_write_->write("critical exponents",exponents[eta_idx]);
			}break;
		default:{ std::cerr<<"std::string AnalyseChain::extract_level_5() : too many wavefunctions"<<std::endl; }
	}
	return filename_;
}

/*different boundary condition*/
std::string AnalyseChain::extract_level_4(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	(*read_)>>nof_;
	/*!must save now nof_ because it doesn't refer to the number of file in the
	 * next directory but in the next-next directory*/
	jd_write_->write("number of different filling",nof_); 

	double polymerization_strength;
	Vector<double> exponents;
	Data<double> E;

	for(unsigned int i(0);i<nof_;i++){
		CreateSystem cs(read_);
		cs.init(read_,this);
		(*read_)>>E>>polymerization_strength>>exponents;

		jd_write_->add_header()->nl();
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

		jd_write_->add_header()->nl();
		cs.save();
	}

	delete read_;
	read_ = NULL;

	return filename_;
}

/*different system size*/
std::string AnalyseChain::extract_level_2(){
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);
	Vector<unsigned int> ref(read.read<Vector<unsigned int> >());
	unsigned int N(read.read<unsigned int>());
	unsigned int m(read.read<unsigned int>());
	
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
	gpenergy.create_image(true);

	Gnuplot gpexp(analyse_+path_+dir_,filename_+"-exponents");
	gpexp+="set multiplot";
	gpexp.margin("0.15","0.85","0.95","0.55");
	gpexp+="unset xtics";
	gpexp.range("x","0","0.037");
	gpexp.range("y","1.3","1.82");
	gpexp+="plot '"+filename_+".dat' u (1/$4):($3==1?$11:1/0) lc 1 t '$\\nu$',\\";
	gpexp+="     "+tostring(2.0-2.0/N)+" lc 1 notitle";
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
	gpexp.create_image(true);

	Gnuplot gppolym(analyse_+path_+dir_,filename_+"-polymerization");
	gppolym.label("x","$n^{-1}$");
	gppolym.label("y2","d-merization strength");
	gppolym.range("x",0,"");
	gppolym.range("y",0,"");
	gppolym+="set key bottom";
	gppolym+="plot '"+filename_+".dat' u (1/$4):($3==1?$9:1/0) t 'P',\\";
	gppolym+="     '"+filename_+".dat' u (1/$4):($3==0?$9:1/0) t 'O'";
	gppolym.save_file();
	gppolym.create_image(true);

	return filename_;
}
