/*!  @file mc.cpp */

#include "ParallelMonteCarlo.hpp"

std::string init(CreateSystem const& cs);
void run(CreateSystem const& cs, double param, std::string path, unsigned int nruns, unsigned int tmax, unsigned int type);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int vec(P.get<unsigned int>("vec"));
	unsigned int type(P.get<unsigned int>("type"));
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	CreateSystem cs(P);
	if(!P.status()){
		std::string path(init(cs));
		switch(vec){
			case 1:{ 
					   Vector<double> param(P.get<Vector<double> >("param"));
					   for(unsigned int i(0);i<param.size();i++){
						   run(cs,param(i),path,nruns,tmax,type);
					   }
				   }break;
			case 2:{ 
					   double param(0);
					   P.get("param",param);
					   run(cs,param,path,nruns,tmax,type);
					   CreateSystem CS(cs,param);
				   }break;
			default:{std::cout<<"unkown simulation"<<std::endl;}
		}
	}
}

std::string init(CreateSystem const& cs){
	std::string path("N"+tostring(cs.get_N())+"/m"+tostring(cs.get_m())+"/n"+tostring(cs.get_n())+"/");
	switch(cs.get_bc()){
		case -1:{path += "A/";}break;
		case 0: {path += "O/";}break;
		case 1: {path += "P/";}break;
		default:{std::cerr<<"Unknown boundary condition"<<std::endl;}
	}
	Linux command;
	command("mkdir -p " + path);
	return path;
}

void run(CreateSystem const& cs, double param, std::string path, unsigned int nruns, unsigned int tmax, unsigned int type){
	CreateSystem CS(cs,param);
	if( CS.use_complex() ){
		ParallelMonteCarlo<std::complex<double> > sim(&CS,path,nruns,tmax,1e9,type);
		sim.run();
	} else {
		ParallelMonteCarlo<double> sim(&CS,path,nruns,tmax,1e9,type);
		sim.run();
	}
}
