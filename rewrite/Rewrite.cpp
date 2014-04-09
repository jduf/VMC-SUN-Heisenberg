#include "Rewrite.hpp"

Rewrite::Rewrite(std::string filename):
	bfile(NULL),
	unlocked(true),
	reading_point_(0)
{
	bfile = fopen(filename.c_str(),"r+");
	if(bfile==NULL){
		unlocked = false;
		std::cerr<<"Rewrite : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

Rewrite::~Rewrite(){
	if(unlocked){
		fclose(bfile);
	} else {
		std::cerr<<"Rewrite : the file "<< filename<< " is locked"<<std::endl;
	}
}

void Rewrite::rewrite_header(std::string s){
	if(unlocked){
		fseek(bfile,0L,SEEK_END);
		unsigned int fsize(ftell(bfile));
		unsigned int N_old(0);
		fseek(bfile,-sizeof(N_old),SEEK_END);
		reading_point_ = fread(&N_old,sizeof(N_old),1,bfile);
		rewind(bfile);
		if(ftruncate(fileno(bfile),fsize-(N_old+sizeof(N_old)))){
			unsigned int N_new(s.size()-1);
			fseek(bfile,0L,SEEK_END);
			fwrite(s.c_str(),1,N_new,bfile);
			fwrite(&N_new,sizeof(N_new),1,bfile);
			fflush(bfile);
		} else {
			std::cerr<<"Rewrite::rewrite_header : the truncation went wrong"<<std::endl;
		}
	} else {
		std::cerr<<"Rewrite : the file "<< filename<< " is locked"<<std::endl;
	}
}
