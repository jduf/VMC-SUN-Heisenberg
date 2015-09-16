#include "Analyse.hpp"

Analyse::Analyse(std::string const& path, unsigned int const& max_level):
	rel_level_(""),
	max_level_(max_level),
	level_(0)
{
	if(path == ""){ study_=0; }
	if(path == "README"){ study_=1; }
	if(path != "README" && path != ""){ study_=2; path_ = path; }
	if(path.find(".jdbin")!=std::string::npos){ study_=3; }
}

void Analyse::do_analyse(){
	switch(study_){
		case 0: /*treat everything*/
			{
				std::cout<<"analysing the whole "<<sim_<<" directory"<<std::endl;
				rst_file_.add_end(std::make_shared<RSTFile>("./","README"));
				IOFiles r("README",false);
				std::string h;
				r>>h;
				rst_file_.first().text(h);
				recursive_search();
				rst_file_.first().save(false,false);
			}break; 
		case 1: /*update only the README file*/
			{
				std::cout<<"analysing only README file"<<std::endl;
				RSTFile rst("./","README");
				IOFiles r("README",false);
				std::string h;
				r>>h;
				rst.text(h);
				Directory d;
				d.search_file_ext(".jdbin",sim_+dir_,false,false);
				d.sort();
				for(unsigned int j(0);j<d.size();j++){
					rst.hyperlink(d.get_name(j),info_+d.get_name(j)+".html");
				}
				rst.save(false,true);
				std::cout<<std::endl<<rst.get()<<std::endl;
			}break;
		case 2: /*treat the directory given as argument*/
			{
				std::cout<<"analysing only directories below "<<path_<<std::endl;
				if(path_[path_.size()-1] != '/'){ path_ += '/'; }
				std::vector<std::string> tmp(my::string_split(path_,'/'));
				path_ = "";
				level_+=1;
				for(unsigned int i(1);i<tmp.size()-1;i++){ 
					path_ += tmp[i]+'/'; 
					level_++;
					rel_level_ += "../";
				}
				dir_ = tmp[tmp.size()-1]+'/';

				rst_file_.add_end(std::make_shared<RSTFile>(info_+path_,dir_.substr(0,dir_.size()-1)));
				recursive_search();
			}break;
		case 3:
			{ std::cerr<<__PRETTY_FUNCTION__<<" : can't analyse a *.jdbin file"<<std::endl; }break;
	}
}

void Analyse::recursive_search(){
	Directory d;
	d.list_dir(sim_+path_+dir_);
	if(d.size()>0){ d.sort(); }
	level_++;
	if(level_>1){ rel_level_ += "../"; }
	for(unsigned int i(0);i<d.size();i++){
		rst_file_.add_end(std::make_shared<RSTFile>(info_+path_+dir_,d.get_name(i)));

		std::string tmp_path(path_);
		std::string tmp_dir(dir_);
		path_ += dir_;
		dir_ = d.get_name(i) + "/";

		if(level_<max_level_){ recursive_search(); }

		path_ = tmp_path;
		dir_ = tmp_dir;
		rst_file_.pop_end();
	}
	search_jdbin();
	level_--;
	rel_level_.erase(0,3);
}

void Analyse::search_jdbin(){
	Directory d;
	d.search_file_ext(".jdbin",sim_+path_+dir_,false,false);
	nof_ = d.size();
	if(d.size()>0){ 
		d.sort();

		Linux command;
		command.mkdir(info_+path_+dir_);
		command.mkdir(analyse_+path_+dir_);
		open_files();

		if(level_==9){
			std::cout<<std::string(6+path_.size()+dir_.size(),' ')<<"|->"<<d.get_name(d.size()-1)<<std::endl;

			filename_ = d.get_name(d.size()-1);
			all_link_names_.add_end(std::make_shared<std::string>(analyse(level_)));
			all_link_files_.add_end(std::make_shared<std::string>((level_==1?info_:dir_)+filename_+".html"));
		} else {
			std::cout<<"lev "<<level_<<" : "<<path_+dir_<<std::endl;
			for(unsigned int i(0); i<d.size();i++){
				std::cout<<std::string(6+path_.size()+dir_.size(),' ')<<"|->"<<d.get_name(i)<<std::endl;

				filename_ = d.get_name(i);
				all_link_names_.add_end(std::make_shared<std::string>(analyse(level_)));
				all_link_files_.add_end(std::make_shared<std::string>((level_==1?info_:dir_)+filename_+".html"));
			}
		}

		all_link_names_.set_target();
		all_link_files_.set_target();
		while ( all_link_names_.target_next() && all_link_files_.target_next() ) {
			rst_file_.last().hyperlink(all_link_names_.get(),all_link_files_.get()); 
		}

		close_files();
		all_link_names_.set();
		all_link_files_.set();
		rst_file_.last().save(false,true);
	}
}

std::string Analyse::extract_level_7(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false);

	System s(*read_);
	CreateSystem cs(&s);
	Vector<double> tmp(read_->read<Vector<double> >());
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);
	/*Only one call of cs.save() is needed*/
	if(!all_link_names_.size()){ 
		s.save_input(*jd_write_);
		jd_write_->add_header()->nl();
		jd_write_->write("number of jdfiles",nof_);
		jd_write_->add_header()->title("System's parameters",'-');
	}
	std::string link_name(cs.analyse(level_));

	delete read_;
	read_ = NULL;

	return link_name;
}
