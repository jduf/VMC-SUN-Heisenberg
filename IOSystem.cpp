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

void IOSystem::set_IOSystem(IOSystem const* const t){ 
	sim_ = t->sim_;
	info_ = t->info_;
	analyse_ = t->analyse_;
	path_ = t->path_;
	dir_ = t->dir_;
	filename_ = t->filename_;
	read_ = t->read_;
	jd_write_ = t->jd_write_;
	data_write_ = t->data_write_;
	rst_file_ = t->rst_file_;
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
		default:
			   {
				   std::cerr<<__PRETTY_FUNCTION__<<" : level="<<level<<" undefined"<<std::endl;
				   return "";
			   }
	}
}
