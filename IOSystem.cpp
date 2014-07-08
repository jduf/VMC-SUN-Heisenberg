#include "IOSystem.hpp"

IOSystem::IOSystem(std::string const& filename):
	sim_("sim/"),
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

IOSystem::IOSystem(IOSystem const& t):
	sim_(t.sim_),
	info_(t.info_),
	analyse_(t.analyse_),
	path_(t.path_),
	dir_(t.dir_),
	filename_(t.filename_),
	read_(t.read_),
	jd_write_(t.jd_write_),
	data_write_(t.data_write_),
	rst_file_(NULL)
{}

void IOSystem::swap_to_assign(IOSystem& t1, IOSystem& t2){
	std::swap(t1.sim_,t2.sim_);
	std::swap(t1.info_,t2.info_);
	std::swap(t1.analyse_,t2.analyse_);
	std::swap(t1.path_,t2.path_);
	std::swap(t1.dir_,t2.dir_);
	std::swap(t1.filename_,t2.filename_);
	std::swap(t1.read_,t2.read_);
	std::swap(t1.jd_write_,t2.jd_write_);
	std::swap(t1.data_write_,t2.data_write_);
}

std::string IOSystem::analyse(unsigned int const& level){
	switch(level){
		case 1:{return extract_level_1();}break;
		case 2:{return extract_level_2();}break;
		case 3:{return extract_level_3();}break;
		case 4:{return extract_level_4();}break;
		case 5:{return extract_level_5();}break;
		case 6:{return extract_level_6();}break;
		default:
			   {
				   std::cerr<<"Analyse::search_jdbin() : level_="<<level<<" undefined"<<std::endl;
				   return "";
			   }
	}
}

void IOSystem::set_IOSystem(IOSystem* t){ *this = *t; }

IOSystem& IOSystem::operator=(IOSystem t){
	swap_to_assign(*this,t);
	return (*this);
}

std::string IOSystem::extract_level_6(){
	return filename_;
}
std::string IOSystem::extract_level_5(){
	return filename_;
}
std::string IOSystem::extract_level_4(){
	return filename_;
}
std::string IOSystem::extract_level_3(){
	return filename_;
}
std::string IOSystem::extract_level_2(){
	return filename_;
}
std::string IOSystem::extract_level_1(){
	return filename_;
}

IOFiles* IOSystem::open_and_get_jd_write(){
	Linux command;
	command("mkdir -p "+path_);
	if(jd_write_){std::cerr<<"already pas bien"<<std::endl;}
	jd_write_ = new IOFiles(path_+filename_+".jdbin",true);
	return jd_write_;
}
