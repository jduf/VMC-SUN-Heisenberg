#include "Read.hpp"

Read::Read(std::string filename, bool binary):
	bfile(NULL),
	locked(false),
	binary(binary)
{
	if(binary){open_binary(filename);}
	else{open_txt(filename);}
}

Read::~Read(){
	if(binary){fclose(bfile);}
	else{tfile.close();}
}

void Read::open_binary(std::string filename){
	filename += ".bin";
	bfile = fopen(filename.c_str(),"rb");
	if(bfile==NULL){
		locked = true;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Read::open_txt(std::string filename){
	filename += ".dat";
	tfile.open(filename.c_str(),std::ios::in);
	if(!tfile){
		locked = true;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	} 
}


