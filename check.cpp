#include "Read.hpp"
#include "Container.hpp"

#include <string>
#include <iostream>
void check(std::string filename);

int main(int argc, char* argv[]){
	if(argc!=2){
		std::cerr<<"check : wrong number of arguements"<<std::endl;
	} else {
		check(argv[1]);
	}
}

void check(std::string filename){
	Read r(filename);
	std::cout<<r.get_header()<<std::endl;
	FileParser file(filename);
	Container param(false);
	std::string wf;
	file.extract<std::string>(wf);
	file.extract<unsigned int>("N",param);
	file.extract<unsigned int>("m",param);
	file.extract<double>("bc",param);

	std::cout<<"wf="<<wf
	<<" N="<<param.get<unsigned int>("N")
	<<" m="<<param.get<unsigned int>("m")
	<<" bc="<<param.get<double>("bc");
	
	if( wf != "chain" ){
		file.extract<unsigned int>("Lx",param);
		file.extract<unsigned int>("Ly",param);
		std::cout<<" Lx="<<param.get<unsigned int>("Lx")
			<<" Ly="<<param.get<unsigned int>("Ly");
		if(wf == "mu" || wf == "trianglemu"){
			file.extract<double>("mu",param);
			std::cout<<" mu="<<param.get<double>("mu");
		}
		if(wf == "phi" || wf == "trianglephi"){
			file.extract<double>("phi",param);
			std::cout<<" phi="<<param.get<double>("phi");
		}
		if( wf == "jastrow" || wf == "trianglejastrow"){
			Vector<double> nu;
			file.extract<Matrix<unsigned int> >("nn",param);
			file.extract<Matrix<unsigned int> >("cc",param);
			file.extract<Vector<unsigned int> >("sl",param);
			file.extract<Matrix<std::complex<double> > >("omega",param);
		}
	}

	file.extract<Matrix<unsigned int> >("sts",param);
	std::cout<<std::endl<<"sts=" <<std::endl
		<<param.get<Matrix<unsigned int> >("sts")<<std::endl;
	if( wf == "mu" || wf == "chain" || wf == "fermi"){
		file.extract<Matrix<double> >("EVec",param);
		std::cout<<"EVec="<<std::endl
			<<param.get<Matrix<double> >("EVec")<<std::endl;
	} 
	if(  wf == "csl" || wf == "phi" || wf == "trianglephi" ){
		file.extract<Matrix<std::complex<double> > >("EVec",param);
		std::cout<<"EVec="<<std::endl
			<<param.get<Matrix<std::complex<double> > >("EVec")<<std::endl;
	}
	if( wf == "jastrow" || wf == "trianglejastrow"){
		std::cout<<"nn="<<std::endl;
		std::cout<<param.get<Matrix<unsigned int> >("nn")<<std::endl;
		std::cout<<"cc="<<std::endl;
		std::cout<<param.get<Matrix<unsigned int> >("cc")<<std::endl;
		std::cout<<"sl="<<std::endl;
		std::cout<<param.get<Vector<unsigned int> >("sl")<<std::endl;
		std::cout<<"omega="<<std::endl;
		std::cout<<param.get<Matrix<std::complex<double> > >("omega")<<std::endl;
	}
}
