#include "Read.hpp"

Read::Read(std::string filename, Matrice<double>& m){
	std::ifstream s;
	s.open(filename.c_str(),std::ios::in);
	if(s){
		for(unsigned int i(0);i<m.size();i++){
			for(unsigned int j(0);j<m.size();j++){
				s>>m(i,j);
			}
		}
	} else{
		std::cerr<<"the file "<< filename<<" doesn't exist"<<std::endl;
	}
	s.close();
}
