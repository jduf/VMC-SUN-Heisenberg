#include "Analyse.hpp"

Analyse::Analyse(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file):
	IOSystem("",sim,"info-"+sim,"analyse-"+sim,"./","./",NULL),
	max_level_(max_level),
	bash_file_(bash_file)
{
	if(sim_ == "README/"){ study_ = 0; }
	else if( std::count(sim_.begin(),sim_.end(),'/') == 1) { study_ = 1; }
	else {
		study_ = 2;

		std::vector<std::string> tmp(my::string_split(sim_,'/'));
		level_=1;
		path_ = "";
		sim_ = tmp[0]+'/';
		info_ = "info-"+sim_;
		analyse_ = "analyse-"+sim_;
		for(unsigned int i(1);i<tmp.size()-1;i++){
			path_ += tmp[i]+'/';
			level_++;
			rel_level_ += "../";
		}
		dir_ = tmp[tmp.size()-1]+'/';
	}
}

Analyse::~Analyse(){ Linux::close(bash_file_>1); }

void Analyse::do_analyse(){
	switch(study_){
		case 0: /*update only the README file*/
			{
				std::cout<<"analysing only README file"<<std::endl;
				RSTFile rst("./","README");
				IOFiles r("README",false,false);
				rst.text(r.read<std::string>());
				rst.save(true,false,true);
			}break;
		case 1: /*treat everything*/
			{
				std::cout<<"analysing the whole "<<sim_<<" directory"<<std::endl;
				list_rst_.add_end(std::make_shared<RSTFile>("./","README"));
				IOFiles r("README",false,false);
				std::string h;
				r>>h;
				list_rst_.first().text(h);
				recursive_search();
				list_rst_.first().save(true,false,false);
			}break;
		case 2: /*only study a given directory and its subdirectories*/
			{
				std::cout<<"analysing only the directories below "<<path_<<std::endl;
				list_rst_.add_end(std::make_shared<RSTFile>(info_+path_,dir_.substr(0,dir_.size()-1)));

				if(level_ != 1 && bash_file_){ Linux::open(dir_.substr(0,dir_.size()-1)+".bash"); }
				recursive_search();
			}break;
	}
}

void Analyse::recursive_search(){
	if(level_ == 1 && bash_file_){ Linux::open(dir_.substr(0,dir_.size()-1)+".bash"); }

	Directory d;
	d.list_dir(sim_+path_+dir_);
	if(d.size()>0){ d.sort(); }
	level_++;
	if(level_>1){ rel_level_ += "../"; }
	for(unsigned int i(0);i<d.size();i++){
		list_rst_.add_end(std::make_shared<RSTFile>(info_+path_+dir_,d.get_name(i)));

		std::string tmp_path(path_);
		std::string tmp_dir(dir_);
		if(dir_ != "./"){ path_ += dir_; }
		dir_ = d.get_name(i) + "/";

		if(level_<max_level_){ recursive_search(); }

		path_ = tmp_path;
		dir_ = tmp_dir;
		list_rst_.pop_end();
	}
	search_jdbin();
	level_--;
	rel_level_.erase(0,3);
}

void Analyse::search_jdbin(){
	Directory d;
	d.search_files_ext(".jdbin",sim_+path_+dir_,false,false);
	nof_ = d.size();
	if(d.size()>0){
		d.sort();

		Linux command;
		command.mkpath((info_+path_+dir_).c_str());
		command.mkpath((analyse_+path_+dir_).c_str());
		open_files();

		if(level_==9 && consider_only_most_recent_jdbin_){
			std::cout<<"lev "<<level_<<" : "<<path_+dir_<<std::endl;
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
		while( all_link_names_.target_next() && all_link_files_.target_next() ){
			list_rst_.last().hyperlink(all_link_names_.get(),all_link_files_.get());
		}

		close_files();
		all_link_names_.set();
		all_link_files_.set();
		list_rst_.last().save(true,false,true);
	}
}

std::string Analyse::extract_default(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	Vector<double> tmp(*read_);
	System s(*read_);

	CreateSystem cs(&s);
	cs.init(&tmp,NULL);
	cs.set_IOSystem(this);

	jd_write_->add_to_header()->nl();
	cs.save(*jd_write_);

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string Analyse::extract_best_of_previous_level(){
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

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
	read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);

	(*read_)>>nof_;
	for(unsigned int i(0);i<idx;i++){
		(*read_)>>param;
		System s(*read_);
	}

	(*read_)>>param;
	System s(*read_);
	CreateSystem cs(&s);
	cs.init(&param,NULL);
	cs.set_IOSystem(this);
	cs.save(*jd_write_);

	delete read_;
	read_ = NULL;

	return filename_;
}

std::string Analyse::fit_therm_limit(){
	if(level_ != 2){ std::cerr<<__PRETTY_FUNCTION__<<" : fit thermodynamical limit is usually done at level 2"<<std::endl; }
	else {
		read_ = new IOFiles(sim_+path_+dir_+filename_+".jdbin",false,false);
		(*read_)>>nof_;

		data_write_ = new IOFiles(analyse_+dir_.substr(0,dir_.size()-1)+".dat",true,false);
		unsigned int ref_check(0);
		for(unsigned int i(0);i<nof_;i++){
			Vector<double> tmp(*read_);
			System s(*read_);
			if(ref_check != s.get_ref()(0)){
				ref_check = s.get_ref()(0);
				std::cout<<"n"<<s.get_n()<<" -> "<<s.get_ref()<<std::endl;
			}
			CreateSystem cs(&s);
			cs.init(&tmp,NULL);
			cs.set_IOSystem(this);
			cs.save(*data_write_);
			if(i+1==nof_){ filename_ = cs.analyse(level_); }
		}

		delete data_write_;
		data_write_ = NULL;

		delete read_;
		read_ = NULL;
	}

	return filename_;
}
