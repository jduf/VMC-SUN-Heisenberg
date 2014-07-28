#include "Analyse.hpp"

Analyse::Analyse(std::string const& sim):
	IOSystem("",sim),
	level_(0)
{}

void Analyse::go(std::string argv){
	Linux command;
	std::string root(command.pwd());
	path_ = argv;
	sim_ = root+sim_;
	info_ = root+info_;
	analyse_ = root+analyse_;
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
				rst_file_.append(RSTFile(root,"README"));
				IOFiles r_readme("README",false);
				std::string h;
				r_readme>>h;
				rst_file_.first().text(h);
				recursive_search();
				rst_file_.first().save(false);
			}break; 
		case 1: /*treat the repository given as argument*/
			{
				if(argv[argv.size()-1] != '/'){ argv += "/"; }
				std::vector<std::string> tmp(string_split(argv,'/'));
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
				analyse_ += tmp[tmp.size()-1] + "/";
				path_ = argv;

				rst_file_.append(RSTFile(info_,filename_));
				recursive_search();
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
	command("mkdir -p " + analyse_+path_+dir_);
	level_++;
	for(unsigned int i(0);i<d.size();i++){
		rst_file_.append(RSTFile(info_+path_+dir_,d.get_name(i)));

		std::string tmp_path(path_);
		std::string tmp_dir(dir_);
		path_ += dir_;
		dir_ = d.get_name(i) + "/";

		recursive_search();

		path_ = tmp_path;
		dir_ = tmp_dir;
		rst_file_.pop();
	}
	search_jdbin();
	level_--;
}

void Analyse::search_jdbin(){
	Directory d;
	d.search_file_ext(".jdbin",sim_+path_+dir_,false,false);
	nof_ = d.size();
	if(d.size()>0){ 
		open_files();

		d.sort();
		std::cout<<"lev "<<level_<<" : "<<sim_+path_+dir_<<std::endl;
		for(unsigned int i(0); i<d.size();i++){
			std::cout<<"-------> "<<d.get_name(i)<<std::endl;
			filename_ = d.get_name(i);
			all_link_names_.append(analyse(level_));
			all_link_files_.append(info_+path_+dir_+filename_+".html");
		}

		for(unsigned int i(0);i<all_link_names_.size();i++){
			rst_file_.last().hyperlink(all_link_names_[i],all_link_files_[i]);
		}

		close_files();
		all_link_names_.clear();
		all_link_files_.clear();
		rst_file_.last().save(false);
	}
}
