#include "Analyse.hpp"

Analyse::Analyse(std::string argv):
	sim_("sim/"),
	info_("info_new/"),
	analysis_("analysis_new/"),
	path_(argv),
	dir_(""),
	filename_(""),
	write_(NULL),
	level_(0)
	//GSR(NULL),
	//GSC(NULL)
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
				RSTFile rst(root,"README");
				IOFiles r_readme("README",false);
				std::string h;
				r_readme>>h;
				rst.text(h);
				recursive_search(rst);
			}break; 
		//case 1: /*treat the repository given as argument*/
			//{
				//if(search_in[search_in.size()-1] != '/'){ search_in += "/"; }
				//info_dir += "info/";
				//std::vector<std::string> tmp(string_split(search_in,'/'));
				//if(tmp.size()<2){
					//std::cerr<<"study : if the update of the whole sim/ directory is requested, then call '\\study' with no argument"<<std::endl;
				//} else {
					//for(unsigned int i(1);i<tmp.size()-1;i++){
						//if(i+2==tmp.size()){/*to update the previous rst file*/
							//Directory d;
							//std::string tmp_local(root+"sim/");
							//for(unsigned int j(1);j<=i;j++){ tmp_local += tmp[j] + "/"; }
							//d.list_dir(tmp_local);
							//d.sort();
							//RSTFile rst(info_dir,tmp[i]);
							//for(unsigned int j(0);j<d.size();j++){
								//rst.hyperlink(d.get_path(j)+d.get_name(j),info_dir+tmp[i]+"/"+d.get_name(j)+".html");
								//rst.nl();
							//}
						//}
						//info_dir += tmp[i] + "/";
						//analysis_dir += tmp[i] + "/";
					//}
					//if(tmp.size()==2){/*to update the previous REAME.rst file*/
						//Directory d;
						//d.list_dir(root+"sim/");
						//d.sort();
						//info_name = "README";
						//IOFiles r_readme(root+"README",false);
						//std::string h;
						//r_readme>>h;
						//RSTFile rst(root,info_name);
						//rst.text(h);
						//for(unsigned int j(0);j<d.size();j++){
							//rst.hyperlink(d.get_path(j)+d.get_name(j),info_dir+d.get_name(j)+".html");
							//rst.nl();
						//}
					//}
					//info_name = tmp[tmp.size()-1];
					//analysis_dir += tmp[tmp.size()-1] + "/";
					//search_in = root + search_in;
//
					//RSTFile rst_next(info_dir,info_name);
					//build_rst(rst_next,search_in,info_dir,analysis_dir,info_name);
				//}
			//}break;
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

void Analyse::print(){
	std::cout<<"------------------"<<std::endl;
	std::cout<<sim_<<std::endl;
	std::cout<<info_<<std::endl;
	std::cout<<analysis_<<std::endl;
	std::cout<<sim_<<std::endl;
	std::cout<<path_<<std::endl;
	std::cout<<dir_<<std::endl;
	std::cout<<"------------------"<<std::endl;
}

void Analyse::recursive_search(RSTFile& rst){
	Directory d;
	d.list_dir(sim_+path_+dir_);
	if(d.size()>0){ d.sort(); }
	Linux command;
	command("mkdir -p " + info_+path_+dir_);
	command("mkdir -p " + analysis_+path_+dir_);
	level_++;
	for(unsigned int i(0);i<d.size();i++){
		RSTFile rst_next(info_+path_+dir_,d.get_name(i));

		std::string tmp_path(path_);
		std::string tmp_dir(dir_);
		path_ += dir_;
		dir_ = d.get_name(i) + "/";
		
		recursive_search(rst_next);

		path_ = tmp_path;
		dir_ = tmp_dir;
	}
	search_jdbin(rst);
	level_--;
}

void Analyse::search_jdbin(RSTFile& rst){
	Directory d;
	d.search_file_ext(".jdbin",sim_+path_+dir_,false,false);
	if(d.size()>0 ){ 
		if(level_>1){ 
			write_ = new IOFiles(sim_+path_+dir_.substr(0,dir_.size()-1)+".jdbin",true);
			unsigned int N(d.size());
			(*write_)("number of jdfiles",N);
			write_->add_to_header("\n");
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
			rst.hyperlink(all_link_names_[i],all_link_files_[i]);
		}

		if(level_ == 5){
			rst.link_figure(analysis_+path_+"P.png","E.png",analysis_+path_+"P.gp",1000);
		}
		if(level_>1){ 
			rst.text(write_->get_header());
			delete write_;
			write_ = NULL;
		}
		all_link_names_.clear();
		all_link_files_.clear();
	}
}

void Analyse::extract_jdbin(){
	switch(level_){
		case 1:{extract_level_1();}break;
		case 2:{extract_level_2();}break;
		case 3:{extract_level_3();}break;
		case 4:{extract_level_4();}break;
		case 5:{extract_sim();}break;
		default:{std::cerr<<"Analyse::search_jdbin(RSTFile& rst) : level_="<<level_<<" undefined"<<std::endl;}
	}
	all_link_files_.append(info_+path_+dir_+filename_+".html");
}

void Analyse::extract_sim(){
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);
	RSTFile rst(info_+path_+dir_,filename_);

	unsigned int type;
	Vector<unsigned int> ref;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	double param;

	read>>type>>ref>>N>>m>>n>>bc>>param;
	CreateSystem CS(ref,N,n,m,bc);
	(*write_)("ref (type of wavefunction)",ref);
	(*write_)("N (N of SU(N))",N);
	(*write_)("m (# of particles per site)",m);
	(*write_)("n (# of site)",n);
	(*write_)("bc (boundary condition)",bc);
	(*write_)("param",param);
	/* the +1 is the averages over all runs */
	CS.treat_one_sim(read,*write_,rst,analysis_+path_+dir_,filename_); 
	rst.text(read.get_header());
	write_->add_to_header("\n");
	all_link_names_.append(tostring(param));
}

void Analyse::extract_level_4(){
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
	IOFiles write(analysis_+path_+dir_+filename_+".dat",true);

	Data<double> min_E;
	min_E.set_x(0.0);
	double min_param(0.0);
	double min_polymerization_strength(0.0);

	for(unsigned int i(0);i<nof;i++){
		read>>ref>>N>>m>>n>>bc>>param>>E>>polymerization_strength;
		write<<param<<" "<<E.get_x()<<" "<<E.get_dx()<<" "<<E.get_conv()<<IOFiles::endl;
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
	(*write_)("min E",min_E);
	(*write_)("param",min_param);
	(*write_)("polymerization strength",min_polymerization_strength);

	Gnuplot gp(analysis_+path_+dir_,filename_);
	gp+="set xlabel '$\\delta$' offset 0,1";
	gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
	gp+="plot 'P.dat' u 1:2:3 w e t 'fix the title issue'";
	gp.save_file();
	gp.create_image(true);
	all_link_names_.append(filename_);
}

void Analyse::extract_level_3(){
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int nof;
	Vector<unsigned int> ref;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	double param;
	Data<double> E;
	double polymerization_strength;
	read>>nof>>ref>>N>>m>>n>>E>>param>>polymerization_strength;
	(*write_)("ref",ref);
	(*write_)("N",N);
	(*write_)("m",m);
	(*write_)("n",n);
	(*write_)("E",E);
	(*write_)("param",param);
	(*write_)("polymerization strength",polymerization_strength);
	all_link_names_.append(filename_);
}

void Analyse::extract_level_2(){
	all_link_names_.append(filename_);
}

void Analyse::extract_level_1(){
	all_link_names_.append(filename_);
}

/*{*/
/*
   virtual void GenericSystem<Type>::treat_one_sim(IOFiles* read, IOFiles* write, std::string const& path, std::string const& filename)=0;
   void ChainPolymerized::treat_one_sim(IOFiles* read, IOFiles* write, std::string const& path, std::string const& filename){
   IOFiles corr_file(info_dir+filename+"-corr.dat",true);
   IOFiles long_range_corr_file(info_dir+filename+"-long-range-corr.dat",true);//should not be delcared when type!=2
   read>>E_>>corr_>>long_range_corr_;
   for(unsigned int j(0);j<corr_.size();j++){
   corr_file<<j+0.5<<" "<<corr_[j]<<IOFiles::endl;
   }
   Vector<double> poly_e_(N/m,0);
   unsigned int i(0);
   while(i<corr.size()){
   for(unsigned int j(0);j<N/m;j++){
   poly_e(j) += corr_(i);
   i++;
   }
   }
   poly_e /= n*m/N;
   poly_e.sort();

   for(unsigned int j(0);j<long_range_corr_.size();j++){
   long_range_corr_file<<j+1<<" "<<long_range_corr_[j]<<IOFiles::endl;
   }
   long_range_corr_file<<IOFiles::endl<<IOFiles::endl;
   if(write){ (*write)<<E<<poly_e(N/m-1)<<poly_e(N/m-2)<<IOFiles::endl; }
   }

   virtual void GenericSystem<Type>::plot_for_one_Nm()=0;
   void ChainPolymerized::plot_for_one_Nm(){
   double min_E(0);
   IOFiles E(info_dir+"E.dat",true);
   double min_polymerization_strength(0.0);
   double polymerization_strength(0.0);
   double min_param(0.0);
   double param(0.0);
   for(unsigned int i(0);i<delta.size();i++){ 
   read>>param;
   E<<param;
   for(unsigned int j(0);i<nruns+1;i++){ 
   read>>E_>>polymerization_strength;
   E<<E_.get_x()<<" "<<E_.get_dx()<<" "<<E_.get_N()<<" "<<E_.get_conv()<<" "<<IOFiles::endl;
   }
   if(E.get_x()<min_E){
   param_min = param;
   min_polymerization_strength = polymerization_strength;
   }
   }

   Gnuplot gp(analysis_dir,"E");
   gp+="set xlabel '$\\delta$' offset 0,1";
   gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
   gp+="plot 'E.dat' u 5:6:7 w e t 'fix the title issue'";
   gp.save_file();
   gp.create_image(true);
   rst.link_figure(analysis_dir+"E.png",analysis_dir+"E.png",analysis_dir+"E.gp",1000);

   (*write)<<N<<m<<n<<para
   }
   */
/*}*/
/*{*/
/*
   virtual void GenericSystem<Type>::treat_one_sim(IOFiles* read, IOFiles* write, std::string const& path, std::string const& filename)=0;
   void ChainPolymerized::treat_one_sim(IOFiles* read, IOFiles* write, std::string const& path, std::string const& filename){
   IOFiles corr_file(info_dir+filename+"-corr.dat",true);
   IOFiles long_range_corr_file(info_dir+filename+"-long-range-corr.dat",true);//should not be delcared when type!=2
   read>>E_>>corr_>>long_range_corr_;
   for(unsigned int j(0);j<corr_.size();j++){
   corr_file<<j+0.5<<" "<<corr_[j]<<IOFiles::endl;
   }
   Vector<double> poly_e_(N/m,0);
   unsigned int i(0);
   while(i<corr.size()){
   for(unsigned int j(0);j<N/m;j++){
   poly_e(j) += corr_(i);
   i++;
   }
   }
   poly_e /= n*m/N;
   poly_e.sort();

   for(unsigned int j(0);j<long_range_corr_.size();j++){
   long_range_corr_file<<j+1<<" "<<long_range_corr_[j]<<IOFiles::endl;
   }
   long_range_corr_file<<IOFiles::endl<<IOFiles::endl;
   if(write){ (*write)<<E<<poly_e(N/m-1)<<poly_e(N/m-2)<<IOFiles::endl; }
   }

   virtual void GenericSystem<Type>::plot_for_one_Nm()=0;
   void ChainPolymerized::plot_for_one_Nm(){
   double min_E(0);
   IOFiles E(info_dir+"E.dat",true);
   double min_polymerization_strength(0.0);
   double polymerization_strength(0.0);
   double min_param(0.0);
   double param(0.0);
   for(unsigned int i(0);i<delta.size();i++){ 
   read>>param;
   E<<param;
   for(unsigned int j(0);i<nruns+1;i++){ 
   read>>E_>>polymerization_strength;
   E<<E_.get_x()<<" "<<E_.get_dx()<<" "<<E_.get_N()<<" "<<E_.get_conv()<<" "<<IOFiles::endl;
   }
   if(E.get_x()<min_E){
   param_min = param;
   min_polymerization_strength = polymerization_strength;
   }
   }

   Gnuplot gp(analysis_dir,"E");
   gp+="set xlabel '$\\delta$' offset 0,1";
   gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
   gp+="plot 'E.dat' u 5:6:7 w e t 'fix the title issue'";
   gp.save_file();
   gp.create_image(true);
   rst.link_figure(analysis_dir+"E.png",analysis_dir+"E.png",analysis_dir+"E.gp",1000);

   (*write)<<N<<m<<n<<para
   }
   */
/*}*/
