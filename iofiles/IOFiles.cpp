#include "IOFiles.hpp"

/*constructors and destructor*/
/*{*/
IOFiles::IOFiles(std::string filename, bool write):
	filename_(filename),
	write_(write),
	binary_(false),
	open_(false),
	file_(NULL),
	header_(NULL)
{
	test_ext();
	if(binary_){open_binary();}
	else{open_txt();}
	if(binary_ && open_){
		if(write_){header_->init(filename_);}
		else{read_header();}
	}
}

IOFiles::~IOFiles(){
	if(open_ && binary_){
		if(write_){write_header();}
		delete header_;
	}
	file_.close();
}

std::string const IOFiles::endl("\n");
/*}*/

/*private methods used in the constructors, destructor or with open(std::string filename_)*/
/*{*/
void IOFiles::test_ext(){
	std::string ext("bin");
	if(filename_.find(ext, (filename_.size() - ext.size())) != std::string::npos){ 
		ext = ".jd" + ext;
		if(filename_.find(ext, (filename_.size() - ext.size())) != std::string::npos){header_ = new Header;}
		binary_ = true;
	}
}

void IOFiles::open_binary(){
	if(write_){file_.open(filename_.c_str(),std::ios::out | std::ios::binary);}
	else {file_.open(filename_.c_str(),std::ios::in | std::ios::binary);}
	if(file_.is_open()){open_ = true;}
	else {std::cerr<<"IOFiles::open_binary() : failed to open "<< filename_<<std::endl;}
}

void IOFiles::open_txt(){
	if(write_){file_.open(filename_.c_str(),std::ios::out);}
	else {file_.open(filename_.c_str(),std::ios::in);}
	if(file_.is_open()){open_ = true;}
	else{std::cerr<<"IOFiles::open_txt() : failed to open "<< filename_<<std::endl;} 
}

void IOFiles::read_header(){
	unsigned int N(0);
	file_.seekg(-sizeof(unsigned int),std::ios::end);
	file_.read((char*)(&N),sizeof(unsigned int));
	char* h(new char[N+1]);
	file_.seekg(-sizeof(char)*N-sizeof(unsigned int),std::ios::end);
	file_.read(h,N);
	h[N] = '\0';
	header_->set(h);
	delete[] h;
	file_.seekg(0,std::ios::beg);
}

void IOFiles::write_header(){
	std::string t(header_->get());
	unsigned int N(t.size());
	file_.write(t.c_str(),N);
	file_.write((char*)(&N),sizeof(unsigned int));
}

void IOFiles::read_string(std::string& t){
	if(open_ && !write_){
		if (binary_){
			unsigned int N(0);
			file_.write((char*)(&N),sizeof(unsigned int));
			char* tmp(new char[N+1]);
			file_.read(tmp,N);
			tmp[N] = '\0';
			t = tmp;
			delete[] tmp;
		} else {
			std::string tmp("");
			while(file_.good()){ 
				getline(file_,tmp);
				t += tmp + "\n";
			}
		}
	} else {
		std::cerr<<"IOFiles::read_basic_type(string) : can't read from "<<filename_<<std::endl;
	}
}

void IOFiles::write_string(const char* t, unsigned int const& N){
	if(open_ && write_){
		if (binary_){
			file_.write((char*)(&N),sizeof(unsigned int));
			file_.write(t,N);
		} else { file_<<t<<std::flush; }
	} else {
		std::cerr<<"IOFiles::write_basic_type(string) : can't write in "<<filename_<<std::endl;
	}
}
/*}*/

/*public methods*/
/*{*/
void IOFiles::precision(unsigned int const& N){ 
	if(binary_){ std::cerr<<"IOFiles::setprecision(unsigned int const& N) : has no effect on a binary file"; }
	else{ file_.precision(N); }
}

std::string IOFiles::get_header() const { 
	if(header_ && open_){
		return header_->get();
	} else {
		std::cerr<<"IOFiles::get_header() : can't read from "<<filename_<<std::endl;
		return 0;
	}
}
/*}*/
