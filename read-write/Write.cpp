#include "Write.hpp"

Write::Write(std::string filename):
	filename(filename),
	bfile(NULL),
	unlocked(true),
	binary(is_binary(filename))
{
	if(binary){open_binary();}
	else{open_txt();}
}

Write::Write():
	filename("no-filename-given"),
	bfile(NULL),
	unlocked(false),
	binary(false)
{}

std::string Write::endl="\n";

void Write::open(std::string filename){
	if(!unlocked){
		this->filename = filename;
		this->binary = is_binary(filename);
		if(binary){open_binary();}
		else{open_txt();}
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

Write::~Write(){
	if(binary){fclose(bfile);}
	else{tfile.close();}
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
