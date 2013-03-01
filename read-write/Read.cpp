#include "Read.hpp"

Read::Read(std::string filename, bool header):
	filename(filename),
	bfile(NULL),
	h(NULL),
	unlocked(true),
	binary(is_binary(filename))
{
	if(binary){open_binary();}
	else{open_txt();}
	if(header && unlocked){h = new Header; read_header();}
}

Read::Read():
	filename("no-filename-given"),
	bfile(NULL),
	h(NULL),
	unlocked(false),
	binary(true)
{ }

Read::~Read(){
	if(unlocked) {
		if(binary){fclose(bfile);}
		else{tfile.close();}
		if(h){ delete h; }
	} else {
		std::cerr<<"Read : the file "<< filename << " is already open"<<std::endl;
	}
}

void Read::open(std::string filename, bool header){
	if(!unlocked){
		this->filename = filename;
		this->binary = is_binary(filename);
		if(binary){open_binary();}
		else{open_txt();}
		if(header){h=new Header;}
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
	std::string ext(".jdbin");
	if(f.find(ext, (f.size() - ext.size())) != std::string::npos){ return true;}
	else{ return false;}
}

void Read::read_header(){
	unsigned int N(0);
	fseek(bfile,-sizeof(N),SEEK_END);
	fread(&N,sizeof(N),1,bfile);
	char c[N];
	fseek(bfile,-sizeof(c)-sizeof(N),SEEK_END);
	fread(c,1,N,bfile);
	c[N] = '\0';
	h->set(c);
	rewind(bfile);
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
			std::string tmp("");
			while(tfile.good()){ 
				getline(tfile,tmp);
				s += tmp + "\n";
			}
		}
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}

