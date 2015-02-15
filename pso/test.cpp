#include"FuncPSO.hpp"
#include"Parseur.hpp"

void run(Parseur& P);

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	run(P);
}

void run(Parseur& P){
	FuncPSO s(P.get<unsigned int>("Nbees"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp"),P.get<unsigned int>("maxiter"));

	s.PSO_init();
	s.PSO_run();
	s.result();
}
