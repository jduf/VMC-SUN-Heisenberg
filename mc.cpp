/*!  @file mc.cpp */

#include "ParallelMonteCarlo.hpp"

//std::string init(CreateSystem const& cs);
template<typename Type>
void run(CreateSystem const& cs, std::string path, unsigned int nruns, unsigned int tmax, unsigned int type);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int type(P.get<unsigned int>("type"));
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	CreateSystem cs(P);
	std::string path("bla");
	if(!P.status()){
		//std::string path(init(cs));
		if(P.is_vector("param")){
			Vector<double> param(P.get<Vector<double> >("param"));
			for(unsigned int i(0);i<param.size();i++){
				cs.create(param(i));
				if( cs.use_complex() ){
					run<std::complex<double> >(cs,path,nruns,tmax,type);
				} else {
					run<double>(cs,path,nruns,tmax,type);
				}
			}
		}
		else {
			double param(P.get<double>("param"));
			cs.create(param);
			if( cs.use_complex() ){
				run<std::complex<double> >(cs,path,nruns,tmax,type);
			} else {
				run<double>(cs,path,nruns,tmax,type);
			}
		}
	}
}

//std::string init(CreateSystem const& cs){
//std::string path("N"+tostring(cs.get_N())+"/m"+tostring(cs.get_m())+"/n"+tostring(cs.get_n())+"/");
//std::string path("bla");
//switch(cs.get_bc()){
//case -1:{path += "A/";}break;
//case 0: {path += "O/";}break;
//case 1: {path += "P/";}break;
//default:{std::cerr<<"Unknown boundary condition"<<std::endl;}
//}
//Linux command;
//command("mkdir -p " + path);
//return path;
//}

template<typename Type>
void run(CreateSystem const& cs, std::string path, unsigned int nruns, unsigned int tmax, unsigned int type){

	//IOFiles results(path+CS.get_filename()+".jdbin",true);
	IOFiles results(path+"bla.jdbin",true);
	results("type of simulation",type);
	results("number of simulations runned",nruns);
	RST rst;
	rst.title("Input","-");
	results.add_to_header(rst.get());
	rst.set();
	cs.save(results);
	rst.title("Results","-");

	MCSystem<Type>* S(NULL);
	if(cs.is_bosonic()){ S = new SystemBosonic<Type>(cs.get_system(),cs.get_bosonic<Type>()); }
	else { S = new SystemFermionic<Type>(cs.get_system(),cs.get_fermionic<Type>()); }
	S->set_type(type);
	ParallelMonteCarlo<Type> sim(S,nruns,tmax);
	sim.run(results);
	sim.save(results);
	delete S;
}
