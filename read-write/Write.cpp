#include "Write.hpp"

Write::Write(std::string filename, bool header):
	filename(filename),
	bfile(NULL),
	h(NULL),
	unlocked(true),
	binary(is_binary(filename))
{
	if(binary){open_binary();}
	else{open_txt();}
	if(header && unlocked ){h = new Header; h->init(filename);}
}

Write::Write():
	filename("no-filename-given"),
	bfile(NULL),
	h(NULL),
	unlocked(false),
	binary(false)
{}

Write::~Write(){
	if(unlocked){
		if(h){write_header(); delete h;}
		if(binary){ fclose(bfile);}
		else{tfile.close();}
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
}

std::string Write::endl="\n";

void Write::open(std::string filename, bool header){
	if(!unlocked){
		this->filename = filename;
		this->binary = is_binary(filename);
		if(binary){open_binary();}
		else{open_txt();}
		if(header){h = new Header; h->set(filename);}
		unlocked=true;
	} else {
		std::cerr<<"Write : the file "<< filename << " is already open"<<std::endl;
	}
}

void Write::open_binary(){
	bfile = fopen(filename.c_str(),"wb");
	if(bfile==NULL){
		unlocked = false;
		std::cerr<<"Write : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Write::open_txt(){
	tfile.open(filename.c_str(),std::ios::out);
	if(!tfile){
		unlocked = false;
		std::cerr<<"Write : the opening of file "<< filename<<" is problematic"<<std::endl;
	} 
}

bool Write::is_binary(std::string f){
	std::string ext(".bin");
	if(f.find(ext, (f.size() - ext.size())) != std::string::npos){ return true;}
	else{ return false;}
}

Write& Write::operator<<(std::string const& s){
	if(unlocked){
		if(binary){
			unsigned int N(s.size());
			fwrite(&N,sizeof(N), 1 ,bfile);
			fwrite(s.c_str(),1, N ,bfile);
			fflush(bfile);
		} else {  
			tfile<<s<<std::flush;
		}
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

void Write::write_header(){
	if(binary){
		if(h){
			std::string s((h->get())->get());
			unsigned int N(s.size());
			fwrite(s.c_str(),1, N ,bfile);
			fwrite(&N,sizeof(N), 1 ,bfile);
			fflush(bfile);
		}
	} else {  
		tfile<<(h->get())->get()<<std::flush;
	}
}
