/*!  @file mc.cpp */

#include "ParallelMonteCarlo.hpp"

std::string init(CreateSystem const& cs);
void run(CreateSystem const& cs, double param, std::string path, unsigned int nruns, unsigned int tmax, unsigned int type);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int type(P.get<unsigned int>("lrc"));
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	CreateSystem cs(P);
	if(!P.status()){
		std::string path(init(cs));
		if(P.is_vector("param")){
			Vector<double> param(P.get<Vector<double> >("param"));
			for(unsigned int i(0);i<param.size();i++){
				run(cs,param(i),path,nruns,tmax,type);
			}
		} else {
			double param(P.get<double>("param"));
			run(cs,param,path,nruns,tmax,type);
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

	IOFiles results(path+CS.get_filename()+".jdbin",true);
	results("type of simulation",type);
	results("number of simulations runned",nruns);
	RST rst;
	rst.title("Input","-");
	results.add_to_header(rst.get());
	rst.set();
	CS.save(results);
	rst.title("Results","-");

	if( CS.use_complex() ){
		ParallelMonteCarlo<std::complex<double> > sim(&CS,nruns,tmax,type);
		sim.run(results);
		sim.save(results);
	} else {
		ParallelMonteCarlo<double> sim(&CS,nruns,tmax,type);
		sim.run(results);
		sim.save(results);
	}
}
