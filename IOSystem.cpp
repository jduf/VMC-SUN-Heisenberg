#include "IOSystem.hpp"

IOSystem::IOSystem(std::string const& filename, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, Vector<unsigned int> const& ref):
	filename_(filename)
{
	std::string tmp;
	tmp = "N" + my::tostring(N);
	filename_ += "-"+tmp;
	path_ += tmp+"/";

	tmp = "m" + my::tostring(m);
	filename_ += "-"+tmp;
	path_ += tmp+"/";

	tmp = "n" + my::tostring(n);
	filename_ += "-"+tmp;
	path_ += tmp+"/";

	tmp = "M";
	for(unsigned int i(0);i<M.size();i++){ tmp  += "_" + my::tostring(M(i)); }
	filename_ += "-"+tmp;
	path_ += tmp+"/";

	switch(bc){
		case-1:{ tmp = "A"; }break;
		case 0:{ tmp = "O"; }break;
		case 1:{ tmp = "P"; }break;
	}
	filename_ += "-"+tmp;
	path_ += tmp+"/";

	tmp = "Juniform";
	filename_ += "-"+tmp;
	path_ += tmp+"/";

	tmp = my::tostring(ref(0))+my::tostring(ref(1))+my::tostring(ref(2));
	path_ += tmp+"/";
}

IOSystem::IOSystem(std::string const& filename, std::string const& sim, std::string const& info, std::string const& analyse, std::string const& path, std::string const& dir, RSTFile* const rst_file, bool const& replace_title_with_link_in_rst):
	sim_(my::ensure_trailing_slash(sim)),
	info_(my::ensure_trailing_slash(info)),
	analyse_(my::ensure_trailing_slash(analyse)),
	path_(my::ensure_trailing_slash(path)),
	dir_(my::ensure_trailing_slash(dir)),
	filename_(filename),
	rst_file_(rst_file),
	replace_title_with_link_in_rst_(replace_title_with_link_in_rst)
{}

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
	replace_title_with_link_in_rst_ = ios->replace_title_with_link_in_rst_;
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
