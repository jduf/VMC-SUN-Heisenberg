#include "Write.hpp"

/*constructors and destructor*/
/*{*/
Write::Write(std::string filename):
	filename(filename),
	bfile(NULL),
	h(NULL),
	unlocked(true),
	binary(test_ext(filename))
{
	if(binary){open_binary();}
	else{open_txt();}
	if(h && unlocked ){h->init(filename);}
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
/*}*/

/*operators*/
/*{*/
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

Write& Write::operator<<(Matrix<std::string> const& mat){
	if(unlocked){
		if(binary){
			unsigned int N_row(mat.row());
			unsigned int N_col(mat.col());
			fwrite(&N_row, sizeof(N_row), 1 ,bfile);
			fwrite(&N_col, sizeof(N_col), 1 ,bfile);
			for(unsigned int i(0); i<N_col*N_row; i++){
				(*this)<<(mat.ptr())[i];
			}
		} else {  
			std::cerr<<"Write : << for Matrix<std::string> is not implemented"<<std::endl;
		}
	} else {
		std::cerr<<"Write : the file "<< filename<< " is locked"<<std::endl;
	}
	return (*this);
}
/*}*/

/*private methods used in the constructors, destructor or with open(std::string filename)*/
/*{*/
void Write::open(std::string filename){
	if(!unlocked){
		this->filename = filename;
		this->binary = test_ext(filename);
		if(binary){open_binary();}
		else{open_txt();}
		if(h){h->init(filename);}
		unlocked=true;
	} else {
		std::cerr<<"Write : the file "<< filename << " is already open"<<std::endl;
	}
}

bool Write::test_ext(std::string f){
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

void Write::write_header(){
	if(binary){
		if(h){
			std::string s((h->get()).get());
			unsigned int N(s.size());
			fwrite(s.c_str(),1, N ,bfile);
			fwrite(&N,sizeof(N), 1 ,bfile);
			fflush(bfile);
		}
	} else {  
		tfile<<(h->get()).get()<<std::flush;
	}
}
/*}*/

/*methods that modify the class*/
/*{*/
void Write::set_header(std::string s){
	if(h){ h->add(s); }
}
/*}*/
