#include "Read.hpp"

Read::Read(std::string filename, bool binary):
	filename(filename),
	bfile(NULL),
	unlocked(true),
	binary(binary)
{
	if(binary){open_binary();}
	else{open_txt();}
}

Read::Read():
	filename("no-filename-given"),
	bfile(NULL),
	unlocked(false),
	binary(binary)
{ }

Read::~Read(){
	if(binary){fclose(bfile);}
	else{tfile.close();}
}

void Read::open(std::string filename, bool binary){
	if(!unlocked){
		this->binary = binary;
		this->filename = filename;
		if(binary){open_binary();}
		else{open_txt();}
		unlocked=true;
	} else {
		std::cerr<<"Read : the file "<< filename << " is already open"<<std::endl;
	}
}

void Read::open_binary(){
	filename += ".bin";
	bfile = fopen(filename.c_str(),"rb");
	if(bfile==NULL){
		unlocked = false;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Read::open_txt(){
	filename += ".dat";
	tfile.open(filename.c_str(),std::ios::in);
	if(!tfile){
		unlocked = false;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	} 
}


