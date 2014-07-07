#include "AnalyseParameter.hpp"


void AnalyseParameter::open_files(std::string const& jdfile, std::string const& datafile, Directory const& d){
	if(level_>1){ 
		jd_write_ = new IOFiles(jdfile,true);
		(*jd_write_)("number of jdfiles",d.size());
		jd_write_->add_to_header("\n");
	}
	if(level_==3 || level_==5){
		data_write_ = new IOFiles(datafile,true);
	}
}

void AnalyseParameter::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 5:{ rst_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","E.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
			case 3:{ rst_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","d-merization.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
		}
		rst_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

void AnalyseParameter::extract_level_5(){

	/*E(param)|n=fixé*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);
	CreateSystem cs(read_);
	cs.init();
	cs.analyse(*this);

	delete read_;
	read_ = NULL;
}

void AnalyseParameter::extract_level_4(){
	/*comparison of E(param_optimal)|n=fixé*/
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int nof(0);
	//double param(0.0);
	//double polymerization_strength;
	(*read_)>>nof;
//
	//Data<double> min_E;
	//min_E.set_x(0.0);
	//double min_param(0.0);
	//double min_polymerization_strength(0.0);
//
	//unsigned int idx(0);
	//for(unsigned int i(0);i<nof;i++){
		//(*read_)>>ref_>>N_>>m_>>n_>>bc_>>param>>E_>>polymerization_strength;
		//if(E_.get_x()<min_E.get_x()){ 
			//idx = i;
			//min_E = E_;
			//min_param = param;
			//min_polymerization_strength = polymerization_strength;
		//}
	//}
//
	//(*write_)("ref (type of wavefunction)",ref_);
	//(*write_)("N (N of SU(N))",N_);
	//(*write_)("m (# of particles per site)",m_);
	//(*write_)("n (# of site)",n_);
	//(*write_)("bc (boundary condition)",bc_);
	//(*write_)("param",param);
	//(*write_)("min E",min_E);
	//(*write_)("polymerization strength",min_polymerization_strength);
//
	//Gnuplot gp(analyse_+path_+dir_,filename_);
	//gp+="set xlabel '$\\delta$' offset 0,1";
	//gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
	//if(idx==0){
		//gp+="f(x) = a+b*x**eta";
		//gp+="a="+tostring(min_E.get_x());
		//gp+="b=1";
		//gp+="eta=1";
		//gp+="set fit quiet";
		//gp+="fit f(x) '"+filename_+".dat' u 1:($4==0?$2:1/0):3 via a,b,eta";
		//gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$'";
		//gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
		//gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean',\\";
		//gp+="     f(x) t sprintf('eta %3.4f',eta)";
	//} else {
		//gp+="f(x) = a+b*(x-c)*(x-c)";
		//gp+="a="+tostring(min_E.get_x());
		//gp+="b=1";
		//gp+="c="+tostring(min_param);
		//gp+="set fit quiet";
		//gp+="fit f(x) '"+filename_+".dat' u 1:($4==0?$2:1/0):3 via a,b,c";
		//gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$'";
		//gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
		//gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean',\\";
		//gp+="     f(x) t sprintf('min %3.4f',c)";
	//}
	//gp.save_file();
	//gp.create_image(true);
	all_link_names_.append(filename_);

	delete read_;
	read_ = NULL;
}

void AnalyseParameter::extract_level_3(){
	/*evolution in function of n*/
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int nof;
	Vector<unsigned int> ref;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	double param;
	Data<double> E;
	double polymerization_strength;
	read>>nof;
	for(unsigned int i(0);i<nof;i++){
		read>>ref>>N>>m>>n>>bc>>param>>E>>polymerization_strength;
		(*data_write_)<<n<<" "<<polymerization_strength<<IOFiles::endl;
		(*jd_write_)("ref",ref);
		(*jd_write_)("N",N);
		(*jd_write_)("m",m);
		(*jd_write_)("n",n);
		(*jd_write_)("bc",bc);
		(*jd_write_)("E",E);
		(*jd_write_)("param",param);
		(*jd_write_)("polymerization strength",polymerization_strength);
	}
	all_link_names_.append(filename_);
}

void AnalyseParameter::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$n^{-1}$'";
	gp.xrange(0,"");
	gp.yrange(0,"");
	gp+="set key bottom";
	gp+="plot '"+filename_+".dat' u (1/$1):2 t 'd-merization strength'";
	gp.save_file();
	gp.create_image(true);
	all_link_names_.append(filename_);
}
