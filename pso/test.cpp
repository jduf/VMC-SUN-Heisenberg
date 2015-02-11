#include"FuncPSO.hpp"
#include"Parseur.hpp"

void run(Parseur& P, std::string filename);
//void rerun(Parseur& P, std::string filename);

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	run(P,"pso.jdbin");
	//rerun(P,"pso.jdbin");
}

void run(Parseur& P, std::string filename){
	FuncPSO s(P.get<unsigned int>("Nbees"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp"),P.get<unsigned int>("maxiter"));

	s.PSO_init();
	s.PSO_run();
	s.result();
	s.PSO_save(filename);
}

//void rerun(Parseur& P, std::string filename){
	//FuncPSO s(P.get<unsigned int>("Nbees"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp"),P.get<unsigned int>("maxiter"));
//
	//s.PSO_load(filename);
	//s.PSO_run();
	//s.result();
//}
