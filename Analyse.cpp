#include "Analyse.hpp"

Analyse::Analyse(std::string argv):
	sim_("sim/"),
	info_("info/"),
	analysis_("analysis/"),
	path_(argv),
	dir_(""),
	filename_(""),
	write_(NULL),
	level_(0)
{
	Linux command;
	std::string root(command.pwd());
	sim_ = root+sim_;
	info_ = root+info_;
	analysis_ = root+analysis_;
	/*info_dir must be root as the RSTFiles are saved in the parent directory*/
	unsigned int study;

	if(path_ == ""){ study=0; }
	else {
		if(path_ != "README"){
			if(path_.find(".jdbin") == std::string::npos){ study = 1; }
			else { study = 2;} 
		} else {study = 3;}
	}

	switch(study){
		case 0: /*treat everything*/
			{
				rst_.append(RSTFile(root,"README"));
				IOFiles r_readme("README",false);
				std::string h;
				r_readme>>h;
				rst_.get()->text(h);
				recursive_search();
				rst_.get()->save(false);
			}break; 
		case 1: /*treat the repository given as argument*/
			{
				if(argv[argv.size()-1] != '/'){ argv += "/"; }
				std::vector<std::string> tmp(string_split(argv,'/'));
				if(tmp.size()<2){
					std::cerr<<"study : if the update of the whole sim/ directory is requested, then call '\\study' with no argument"<<std::endl;
				} else {
					for(unsigned int i(1);i<tmp.size()-1;i++){
						if(i+2==tmp.size()){/*to update the previous rst file*/
							Directory d;
							std::string tmp_local(sim_);
							for(unsigned int j(1);j<=i;j++){ tmp_local += tmp[j] + "/"; }
							d.list_dir(tmp_local);
							d.sort();
							RSTFile rst(info_,tmp[i]);
							for(unsigned int j(0);j<d.size();j++){
								rst.hyperlink(d.get_path(j)+d.get_name(j),info_+tmp[i]+"/"+d.get_name(j)+".html");
								rst.nl();
							}
						}
						path_ += tmp[i] + "/";
					}
					if(tmp.size()==2){/*to update the previous REAME.rst file*/
						Directory d;
						d.list_dir(root+"sim/");
						d.sort();
						filename_ = "README";
						IOFiles rst_readme(root+filename_,false);
						std::string h;
						rst_readme>>h;
						RSTFile rst(root,filename_);
						rst.text(h);
						for(unsigned int j(0);j<d.size();j++){
							rst.hyperlink(d.get_path(j)+d.get_name(j),info_+d.get_name(j)+".html");
							rst.nl();
						}
					}
					info_ = tmp[tmp.size()-1];
					analysis_ += tmp[tmp.size()-1] + "/";
					path_ = argv;

					rst_.append(RSTFile(info_,filename_));
					recursive_search();
				}
			}break;
			//case 2:  /*treat only one jdbin file*/
			//{
			//info_dir += "info/";
			//const std::string ext(".jdbin");
			//std::string filename(search_in.substr(0, search_in.size() - ext.size()));
			//std::string path(root);
			//std::vector<std::string> tmp(string_split(filename,'/'));
			//for(unsigned int i(0);i<tmp.size()-1;i++){
			//path += tmp[i] + "/";
			//}
			//for(unsigned int i(1);i<tmp.size()-1;i++){
			//info_dir += tmp[i] + "/";
			//}
			//filename = tmp[tmp.size()-1];
			//
			//extract_jdbin(path,path,filename);
			//}break;
			//case 3: /*update only the README file*/
			//{
			//RSTFile rst(root,"README");
			//IOFiles r("README",false);
			//std::string h;
			//r>>h;
			//rst.text(h);
			//Directory d;
			//d.list_dir(root+"sim/");
			//d.sort();
			//info_dir += "info/";
			//for(unsigned int j(0);j<d.size();j++){
			//rst.hyperlink(d.get_path(j)+d.get_name(j),info_dir+d.get_name(j)+".html");
			//rst.nl();
			//}
			//}break;
	}
}

void Analyse::recursive_search(){
	Directory d;
	d.list_dir(sim_+path_+dir_);
	if(d.size()>0){ d.sort(); }
	Linux command;
	command("mkdir -p " + info_+path_+dir_);
	command("mkdir -p " + analysis_+path_+dir_);
	level_++;
	for(unsigned int i(0);i<d.size();i++){
		rst_.append(RSTFile(info_+path_+dir_,d.get_name(i)));

		std::string tmp_path(path_);
		std::string tmp_dir(dir_);
		path_ += dir_;
		dir_ = d.get_name(i) + "/";

		recursive_search();

		path_ = tmp_path;
		dir_ = tmp_dir;
		rst_.pop();
	}
	search_jdbin();
	level_--;
}

void Analyse::search_jdbin(){
	Directory d;
	d.search_file_ext(".jdbin",sim_+path_+dir_,false,false);
	if(d.size()>0){ 
		if(level_>1){ 
			write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
			(*write_)("number of jdfiles",d.size());
			write_->add_to_header("\n");
		}
		if(level_>4){
			data_write_ = new IOFiles(analysis_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);
		}

		d.sort();
		std::cout<<"lev : "<<level_<<" rep : "<<sim_+path_+dir_<<std::endl;
		for(unsigned int i(0); i<d.size();i++){
			std::cout<<"----->"<<d[i]<<std::endl;
			filename_ = d.get_name(i);
			extract_jdbin();
		}

		//Vector<unsigned int> index(all_link_names_.sort());
		//all_link_files_ = all_link_files_.sort(index);
		for(unsigned int i(0);i<all_link_names_.size();i++){
			rst_.last()->hyperlink(all_link_names_[i],all_link_files_[i]);
		}

		if(level_>1){ 
			if(level_ == 5){
				rst_.last()->link_figure(analysis_+path_+dir_.substr(0,dir_.size()-1)+".png","E.png",analysis_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000);
			}
			rst_.last()->text(write_->get_header());
			delete write_;
			write_ = NULL;
		}
		if(level_>4){
			delete data_write_;
			data_write_ = NULL;
		}
		all_link_names_.clear();
		all_link_files_.clear();
		rst_.last()->save(false);
	}
}

void Analyse::extract_jdbin(){
	switch(level_){
		case 1:{extract_level_1();}break;
		case 2:{extract_level_2();}break;
		case 3:{extract_level_3();}break;
		case 4:{extract_level_4();}break;
		case 5:{extract_sim();}break;
		default:{std::cerr<<"Analyse::search_jdbin() : level_="<<level_<<" undefined"<<std::endl;}
	}
	all_link_files_.append(info_+path_+dir_+filename_+".html");
}

void Analyse::extract_sim(){
	/*E(param)|n=fixé*/
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);
	RSTFile rst(info_+path_+dir_,filename_);

	unsigned int nruns;
	double param;

	read>>nruns>>ref_>>N_>>m_>>n_>>bc_>>param;
	IOFiles corr_file(analysis_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles long_range_corr_file(analysis_+path_+dir_+filename_+"-long-range-corr.dat",true);//should not be delcared when type!=2
	data_write_->precision(10);
	(*data_write_)<<"% delta E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		read>>E_>>corr_>>long_range_corr_;
		(*data_write_)<<param<<" "<<E_.get_x()<<" "<<E_.get_dx();
		if(i<nruns){ (*data_write_)<<" 0"<<IOFiles::endl; }
		else { (*data_write_)<<" 1"<<IOFiles::endl; }
		for(unsigned int j(0);j<corr_.size();j++){
			corr_file<<j+0.5<<" "<<corr_[j]<<IOFiles::endl;
		}
		corr_file<<IOFiles::endl<<IOFiles::endl;
		for(unsigned int j(0);j<long_range_corr_.size();j++){
			long_range_corr_file<<j+1<<" "<<long_range_corr_[j]<<IOFiles::endl;
		}
		long_range_corr_file<<IOFiles::endl<<IOFiles::endl;
	}

	Vector<double> poly_e(N_/m_,0);
	unsigned int i(0);
	while(i<corr_.size()){
		for(unsigned int j(0);j<N_/m_;j++){
			poly_e(j) += corr_[i].get_x();
			i++;
		}
	}
	poly_e /= n_*m_/N_;
	poly_e.sort();

	(*write_)("ref (type of wavefunction)",ref_);
	(*write_)("N (N of SU(N))",N_);
	(*write_)("m (# of particles per site)",m_);
	(*write_)("n (# of site)",n_);
	(*write_)("bc (boundary condition)",bc_);
	(*write_)("param",param);
	(*write_)("E",E_);
	(*write_)("polymerization strength",poly_e(N_/m_-2)-poly_e(N_/m_-1));

	/*{*/
	Gnuplot gp(analysis_+path_+dir_,filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(param)+"$'";
	gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
	gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
	gp+="                   '"+filename_+"-corr.dat' i " + tostring(nruns) + " u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle";
	gp.save_file();
	rst.link_figure(analysis_+path_+dir_+filename_+"-corr.png","Correlation on links",analysis_+path_+dir_+filename_+"-corr.gp",1000);

	unsigned int length(long_range_corr_.size());
	if(length>0){
		Gnuplot gp(analysis_+path_+dir_,filename_+"-long-range-corr");
		gp+="stats '"+filename_+"-long-range-corr.dat' nooutput";
		gp.xrange(0,length+1);
		gp.yrange("1.1*STATS_min_y","1.1*STATS_max_y");
		gp+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(param)+"$'";
		gp+="set key right bottom";
		gp+="a=1.0";
		gp+="b=1.0";
		gp+="eta=1.0";
		gp+="m="+tostring(m_)+".0";
		gp+="N="+tostring(N_)+".0";
		gp+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
		gp+="set fit quiet";
		switch(N_/m_){
			case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; } break;
			case 3:{
					   switch((length + 1) % 3){
						   case 0:{ gp+="fit [2:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 1:{ gp+="fit [5:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   }break;
				   }break;
			default :{ gp+="fit ["+tostring(N_-1)+":] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
		}
		gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-long-range-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
		gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-long-range-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
		gp+="                   '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
		gp+="                   f(x) notitle";
		gp.save_file();
		rst.link_figure(analysis_+path_+dir_+filename_+"-long-range-corr.png","Long range correlation",analysis_+path_+dir_+filename_+"-long-range-corr.gp",1000);
	}
	/*}*/

	rst.text(read.get_header());
	rst.save(false);
	all_link_names_.append(tostring(param));
}

void Analyse::extract_level_4(){
	/*comparison of E(param_optimal)|n=fixé*/
	//IOFiles read(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",false);
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int nof(0);
	Vector<unsigned int> ref;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc(0);
	double param(0.0);
	Data<double> E;
	double polymerization_strength;
	read>>nof;
	//IOFiles write(analysis_+path_+dir_.substr(0,dir_.size()-1)+".dat",true);

	Data<double> min_E;
	min_E.set_x(0.0);
	double min_param(0.0);
	double min_polymerization_strength(0.0);

	for(unsigned int i(0);i<nof;i++){
		read>>ref>>N>>m>>n>>bc>>param>>E>>polymerization_strength;
		if(E.get_x()<min_E.get_x()){ 
			min_E = E;
			min_param = param;
			min_polymerization_strength=polymerization_strength;
		}
	}
	(*write_)("ref",ref);
	(*write_)("N",N);
	(*write_)("m",m);
	(*write_)("n",n);
	(*write_)("bc",bc);
	(*write_)("param",min_param);
	(*write_)("min E",min_E);
	(*write_)("polymerization strength",min_polymerization_strength);

	Gnuplot gp(analysis_+path_+dir_,filename_);
	gp+="set xlabel '$\\delta$' offset 0,1";
	gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
	gp+="f(x) = a+b*(x-c)*(x-c)";
	gp+="a=-1";
	gp+="b=1";
	gp+="c=1";
	gp+="set fit quiet";
	gp+="fit f(x) '"+filename_+".dat' u 1:($4==1?$2:1/0):3 via a,b,c";
	gp+="set title '$N="+tostring(N)+"$ $m="+tostring(n)+"$ $n="+tostring(n)+"$'";
	gp+="plot '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Independant measures',\\";
	gp+="     '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Mean',\\";
	gp+="     f(x) t sprintf('min %3.4f',c)";
	gp.save_file();
	gp.create_image(true);
	all_link_names_.append(filename_);
}

void Analyse::extract_level_3(){
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
		(*write_)("ref",ref);
		(*write_)("N",N);
		(*write_)("m",m);
		(*write_)("n",n);
		(*write_)("bc",bc);
		(*write_)("E",E);
		(*write_)("param",param);
		(*write_)("polymerization strength",polymerization_strength);
	}
	all_link_names_.append(filename_);
}

void Analyse::extract_level_2(){
	all_link_names_.append(filename_);
}

void Analyse::extract_level_1(){
	all_link_names_.append(filename_);
}
