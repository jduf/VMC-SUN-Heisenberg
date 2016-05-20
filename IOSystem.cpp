#include "IOSystem.hpp"

IOSystem::IOSystem(std::string const& filename, std::vector<std::string> names):
	filename_(filename)
{
	for(unsigned int i(0);i<names.size()-1;i++){
		filename_ += "-"+names[i];
		path_ += names[i]+"/";
	}
	path_ += names.back()+"/";
}

IOSystem::IOSystem(std::string const& filename, std::string const& sim, std::string const& info, std::string const& analyse, std::string const& path, std::string const& dir, RSTFile* const rst_file):
	sim_(sim),
	info_(info),
	analyse_(analyse),
	path_(path),
	dir_(dir),
	filename_(filename),
	rst_file_(rst_file)
{
	if(sim_.back()!='/'){ sim_ += "/"; }
	if(info_.back()!='/'){ info_ += "/"; }
	if(analyse_.back()!='/'){ analyse_ += "/"; }
}

void IOSystem::set_IOSystem(IOSystem const* const ios){
	sim_       = ios->sim_;
	info_      = ios->info_;
	analyse_   = ios->analyse_;
	path_      = ios->path_;
	dir_       = ios->dir_;
	filename_  = ios->filename_;
	read_      = ios->read_;
	jd_write_  = ios->jd_write_;
	data_write_= ios->data_write_;
	rst_file_  = ios->rst_file_;
}

std::string IOSystem::analyse(unsigned int const& level){
	switch(level){
		case 1:{ return extract_level_1(); }break;
		case 2:{ return extract_level_2(); }break;
		case 3:{ return extract_level_3(); }break;
		case 4:{ return extract_level_4(); }break;
		case 5:{ return extract_level_5(); }break;
		case 6:{ return extract_level_6(); }break;
		case 7:{ return extract_level_7(); }break;
		case 8:{ return extract_level_8(); }break;
		case 9:{ return extract_level_9(); }break;
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : level="<<level<<" undefined"<<std::endl;
	return "";
}
