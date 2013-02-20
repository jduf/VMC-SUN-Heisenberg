#include "Write.hpp"

Write::Write(std::string filename, bool binary):
	filename(filename),
	bfile(NULL),
	unlocked(true),
	binary(binary)
{
	if(binary){open_binary();}
	else{open_txt();}
}

Write::Write():
	filename("no-filename-given"),
	bfile(NULL),
	unlocked(false),
	binary(binary)
{}

std::string Write::endl="\n";

void Write::open(std::string filename, bool binary){
	if(!unlocked){
		this->binary = binary;
		this->filename = filename;
		if(binary){open_binary();}
		else{open_txt();}
		unlocked=true;
	} else {
		std::cerr<<"Write : the file "<< filename << " is already open"<<std::endl;
	}
}

void Write::open_binary(){
	filename += ".bin";
	bfile = fopen(filename.c_str(),"wb");
	if(bfile==NULL){
		unlocked = false;
		std::cerr<<"Write : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Write::open_txt(){
	filename += ".dat";
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


