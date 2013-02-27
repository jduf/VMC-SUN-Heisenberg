#include "Read.hpp"

Read::Read(std::string filename):
	filename(filename),
	bfile(NULL),
	unlocked(true),
	binary(is_binary(filename))
{
	if(binary){open_binary();}
	else{open_txt();}
}

Read::Read():
	filename("no-filename-given"),
	bfile(NULL),
	unlocked(false),
	binary(true)
{ }

Read::~Read(){
	if(binary){fclose(bfile);}
	else{tfile.close();}
}


void Read::open(std::string filename){
	if(!unlocked){
		this->filename = filename;
		this->binary = is_binary(filename);
		if(binary){open_binary();}
		else{open_txt();}
		unlocked=true;
	} else {
		std::cerr<<"Read : the file "<< filename << " is already open"<<std::endl;
	}
}

void Read::open_binary(){
	bfile = fopen(filename.c_str(),"rb");
	if(bfile==NULL){
		unlocked = false;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	}

}

void Read::open_txt(){
	tfile.open(filename.c_str(),std::ios::in);
	if(!tfile){
		unlocked = false;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	} 
}

bool Read::is_binary(std::string f){
	std::string ext(".bin");
	if(f.find(ext, (f.size() - ext.size())) != std::string::npos){ return true;}
	else{ return false;}
}

Read& Read::operator>>(std::string& s){
	if(unlocked){
		if(binary){ 
			unsigned int N(0);
			fread(&N,sizeof(N),1,bfile);
			char c[N];
			fread(c,1,N,bfile);
			c[N] = '\0';
			s = c;
		} else {
			std::cerr<<"Read : << for string in text file is not implemented"<<std::endl;
		}
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}
