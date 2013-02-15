#include "Write.hpp"

Write::Write(std::string filename, bool binary):
	bfile(NULL),
	binary(binary),
	locked(false)
{
	if(binary){open_binary(filename);}
	else{open_txt(filename);}
}

void Write::open_binary(std::string filename){
	bfile = fopen(filename.c_str(),"wb");
	if(bfile==NULL){
		locked = true;
		std::cerr<<"Write : the opening of file "<< filename<<" is problematic"<<std::endl;
	}
}

void Write::open_txt(std::string filename){
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

Write& Write::operator<<(Matrice<double> const& m){
	if(binary) { write_binary_matrix(m); }
	else { write_txt_matrix(m); }

	return (*this);
}

void Write::write_binary_matrix(Matrice<double> const& m){
	unsigned int N(m.size());
	fwrite(&N, sizeof(N), 1 ,bfile);
	double tmp[N*N];
	for(unsigned int i(0);i<N*N;i++){
		tmp[i] = (m.ptr())[i];
	}
	fwrite(&tmp, sizeof(tmp), 1 ,bfile);
	fflush(bfile);
}

void Write::write_txt_matrix(Matrice<double> const& m){
	for(unsigned int i(0); i<m.size();i++){
		for(unsigned int j(0); j<m.size();j++){
			tfile << m(i,j)<<" ";
		}
		tfile<<std::endl;
	}
}
