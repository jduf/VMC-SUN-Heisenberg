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
	bfile = fopen(filename.c_str(),"rb");
	if(bfile==NULL){
		locked = true;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Read::open_txt(std::string filename){
	tfile.open(filename.c_str(),std::ios::in);
	if(!tfile){
		locked = true;
		std::cerr<<"Read : the opening of file "<< filename<<" is problematic"<<std::endl;
	} 
}

Read& Read::operator>>(Matrice<double>& m){
	if(binary) { read_binary_matrix(m); }
	else { read_txt_matrix(m); }
	return (*this);
}

void Read::read_binary_matrix(Matrice<double>& m){
	unsigned int N(0);
	fread(&N,sizeof(N),1,bfile);
	if(N != m.size()) {
		Matrice<double> mat_tmp(N);
		double tmp[N*N];
		fread(&tmp,sizeof(tmp),1,bfile);
		for(unsigned int i(0);i<N*N;i++){
			(mat_tmp.ptr())[i]=tmp[i];
		}
		m = mat_tmp;
	} else {
		double tmp[N*N];
		fread(&tmp,sizeof(tmp),1,bfile);
		for(unsigned int i(0);i<N*N;i++){
			(m.ptr())[i]=tmp[i];
		}
	}
}

void Read::read_txt_matrix(Matrice<double>& m){
	if(m.size()!=0) {
		for(unsigned int i(0); i<m.size();i++){
			for(unsigned int j(0); j<m.size();j++){
				tfile >> m(i,j);
			}
		}
	} else {
		std::cerr<<"Read : to read a Matrice<double> you need to set the size of the input matric"<<std::endl;
	}
}
