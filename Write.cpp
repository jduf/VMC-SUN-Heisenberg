#include "Write.hpp"

Write::Write(std::string filename, bool binary):
	bfile(NULL),
	binary(binary),
	locked(false)
{
	if(binary){open_binary(filename);}
	else{open_txt(filename);}
}

std::string Write::endl="\n";

void Write::open_binary(std::string filename){
	filename += ".bin";
	bfile = fopen(filename.c_str(),"wb");
	if(bfile==NULL){
		locked = true;
		std::cerr<<"Write : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Write::open_txt(std::string filename){
	filename += ".dat";
	tfile.open(filename.c_str(),std::ios::out);
	if(!tfile){
		locked = true;
		std::cerr<<"Write : the opening of file "<< filename<<" is problematic"<<std::endl;
	} 
}

Write::~Write(){
	if(binary){fclose(bfile);}
	else{tfile.close();}
}


