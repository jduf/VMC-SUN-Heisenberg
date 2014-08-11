#include "IOSystem.hpp"

IOSystem::IOSystem(std::string const& filename, std::string const& sim):
	sim_(sim),
	info_("info/"),
	analyse_("analyse/"),
	path_(""),
	dir_(""),
	filename_(filename),
	read_(NULL),
	jd_write_(NULL),
	data_write_(NULL),
	rst_file_(NULL)
{}

void IOSystem::set_IOSystem(IOSystem* t){ 
	sim_ = t->sim_;
	info_ = t->info_;
	analyse_ = t->analyse_;
	path_ = t->path_;
	dir_ = t->dir_;
	filename_ = t->filename_;
	read_ = t->read_;
	jd_write_ = t->jd_write_;
	data_write_ = t->data_write_;
}

std::string IOSystem::get_filename() const { 
	Linux command;
	command("/bin/mkdir -p "+path_);
	return path_+filename_; 
}

void IOSystem::init_output_file(IOFiles& output){ 
	if(jd_write_){std::cerr<<"void IOSystem::init_output_file(IOFiles& output) : a jd_write_ has alrady been declared"<<std::endl;}
	else{
		jd_write_ = &output; 
		RST rst;
		rst.title("Input","-");
		jd_write_->add_to_header(rst.get());
	}
}

std::string IOSystem::analyse(unsigned int const& level){
	switch(level){
		case 1:{return extract_level_1();}break;
		case 2:{return extract_level_2();}break;
		case 3:{return extract_level_3();}break;
		case 4:{return extract_level_4();}break;
		case 5:{return extract_level_5();}break;
		case 6:{return extract_level_6();}break;
		case 7:{return extract_level_7();}break;
		default:
			   {
				   std::cerr<<"IOSystem::analyse(unsigned int const& level) : level="<<level<<" undefined"<<std::endl;
				   return "";
			   }
	}
}
