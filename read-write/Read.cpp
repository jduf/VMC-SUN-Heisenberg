#include "Read.hpp"

/*constructors and destructor*/
/*{*/
Read::Read(std::string filename):
	filename(filename),
	bfile(NULL),
	h(NULL),
	unlocked(true),
	binary(test_ext(filename)),
	reading_point_(0),
	size_without_header(0)
{
	if(binary){open_binary();}
	else{open_txt();}
	if(h){h->set(read_header());}
}

Read::Read():
	filename("no-filename-given"),
	bfile(NULL),
	h(NULL),
	unlocked(false),
	binary(true),
	reading_point_(0)
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
/*}*/

/*operators*/
/*{*/
Read& Read::operator>>(std::string& s){
	if(unlocked){
		if(binary){ 
			unsigned int N(0);
			reading_point_ = fread(&N,sizeof(N),1,bfile);
			char* c(new char[N+1]);
			reading_point_ = fread(c,1,N,bfile);
			c[N] = '\0';
			s = c;
			delete[] c;
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

Read& Read::operator>>(Matrix<std::string>& mat){
	if(unlocked){
		if(binary){ 
			unsigned int N_row(0);
			unsigned int N_col(0);
			reading_point_ = fread(&N_row,sizeof(N_row),1,bfile);
			reading_point_ = fread(&N_col,sizeof(N_col),1,bfile);
			if(mat.row() != N_row || mat.col() != N_col){ mat.set(N_row,N_col); } 
			for(unsigned int i(0);i<N_row*N_col;i++){
				(*this)>>(mat.ptr())[i];
			}			
		} else {
			std::cerr<<"Read : << for Matrix<std::string> is not implemented"<<std::endl;
		}
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}
/*}*/

/*private methods used in the constructors or with open(std::string filename)*/
/*{*/
void Read::open(std::string filename){
	if(!unlocked){
		this->filename = filename;
		this->binary = test_ext(filename);
		if(binary){open_binary();}
		else{open_txt();}
		unlocked=true;
		if(h){read_header();}
	} else {
		std::cerr<<"Read : the file "<< filename << " is already open"<<std::endl;
	}
}

bool Read::test_ext( std::string f){
	std::string base_ext("bin");
	std::string ext("");
	if(f.find(ext, (f.size() - ext.size())) != std::string::npos){ 
		ext = "." + base_ext;
		if(f.find(ext, (f.size() - ext.size())) != std::string::npos){ return true;}
		ext = ".jd" + base_ext;
		if(f.find(ext, (f.size() - ext.size())) != std::string::npos){ h = new Header; return true; }
	}
	return false;
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

std::string Read::read_header(){
	std::string header("");
	if(unlocked){ 
		if(binary){
			unsigned int N(0);
			fseek(bfile,-sizeof(N),SEEK_END);
			reading_point_ = fread(&N,sizeof(N),1,bfile);
			char* c(new char[N+1]);
			fseek(bfile,-sizeof(char)*N-sizeof(N),SEEK_END);
			size_without_header = ftell(bfile);
			reading_point_ = fread(c,1,N,bfile);
			c[N] = '\0';
			rewind(bfile);
			header = c;
			delete[] c;
		} else {
			std::cerr<<"Read : read_header() not implemented for non binary file"<<std::endl;
		}
	} else {
		std::cerr<<"Read : the file "<< filename<< " is locked"<<std::endl;
	}
	return header;
}
/*}*/

/*methods that returns something related to the class*/
/*{*/
std::string Read::get_header() const { 
	if(h && unlocked){
		return h->get();
	} else {
		std::cerr<<"Read : the file has no header or is locked"<<std::endl;
		return 0;
	}
}

bool Read::get_status() const{
	return unlocked;
}
/*}*/
