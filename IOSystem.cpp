#include "IOSystem.hpp"

IOSystem::IOSystem(std::string const& filename):
	sim_(""),
	info_(""),
	analyse_(""),
	path_(""),
	dir_(""),
	filename_(filename),
	read_(NULL),
	jd_write_(NULL),
	data_write_(NULL)
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
	data_write_(t.data_write_)
{}

void IOSystem::swap_to_assign(IOSystem& t1, IOSystem& t2){
	std::swap(t1.sim_,t2.sim_);
	std::swap(t1.info_,t2.info_);
	std::swap(t1.analyse_,t2.analyse_);
	std::swap(t1.path_,t2.path_);
	std::swap(t1.filename_,t2.filename_);
	std::swap(t1.read_,t2.read_);
	std::swap(t1.jd_write_,t2.jd_write_);
	std::swap(t1.data_write_,t2.data_write_);
}

void IOSystem::analyse(IOSystem const& t){
	*this = t;
}

IOSystem& IOSystem::operator=(IOSystem t){
	swap_to_assign(*this,t);
	return (*this);
}
